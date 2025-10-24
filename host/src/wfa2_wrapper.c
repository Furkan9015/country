/**
 * WFA2 Wavefront Alignment Wrapper for UPVC
 * Replaces O(n²) Smith-Waterman with O(s) WFA2 for 10-20× speedup
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wfa2_wrapper.h"
#include "upvc.h"

// WFA2 library headers
#include "wavefront/wavefront_align.h"

/**
 * Convert UPVC's 2-bit encoded nucleotides to ASCII string
 * UPVC encoding: A=0, C=1, T=2, G=3 (stored in bits 0-1 of int8_t)
 */
static void upvc_to_ascii(const int8_t *encoded, int len, char *ascii) {
    static const char nucleotides[] = "ACTG";
    for (int i = 0; i < len; i++) {
        ascii[i] = nucleotides[encoded[i] & 0x3];
    }
    ascii[len] = '\0';
}

/**
 * Convert WFA2 CIGAR to UPVC backtrack format
 *
 * WFA2 CIGAR operations:
 * - 'M' = Match/Mismatch
 * - 'X' = Mismatch
 * - '=' = Match
 * - 'I' = Insertion (in read relative to reference)
 * - 'D' = Deletion (in read relative to reference)
 *
 * UPVC backtrack:
 * - type: 0=match, 1=substitution, 2=insertion, 3=deletion
 * - ix, jx: positions in reference and read
 */
static int cigar_to_backtrack(cigar_t *cigar, const int8_t *ref, const int8_t *read,
                               backtrack_t *backtrack, int size_neighbour_in_symbols) {
    int backtrack_idx = 0;
    int ref_pos = 0;
    int read_pos = 0;

    // Process CIGAR in reverse (UPVC stores backtrack from end to start)
    for (int i = cigar->begin_offset; i < cigar->end_offset; i++) {
        char op = cigar->operations[i];

        switch (op) {
            case 'M':  // Match or mismatch - need to check
            case '=':  // Explicit match
            case 'X':  // Explicit mismatch
                if ((ref[ref_pos] & 0x3) != (read[read_pos] & 0x3)) {
                    // Substitution
                    backtrack[backtrack_idx].type = 1;  // CODE_SUB equivalent
                    backtrack[backtrack_idx].ix = ref_pos;
                    backtrack[backtrack_idx].jx = read_pos;
                    backtrack_idx++;
                }
                ref_pos++;
                read_pos++;
                break;

            case 'I':  // Insertion in read
                backtrack[backtrack_idx].type = 2;  // CODE_INS equivalent
                backtrack[backtrack_idx].ix = ref_pos;
                backtrack[backtrack_idx].jx = read_pos;
                backtrack_idx++;
                read_pos++;
                break;

            case 'D':  // Deletion in read
                backtrack[backtrack_idx].type = 3;  // CODE_DEL equivalent
                backtrack[backtrack_idx].ix = ref_pos;
                backtrack[backtrack_idx].jx = read_pos;
                backtrack_idx++;
                ref_pos++;
                break;
        }
    }

    // Add end marker
    backtrack[backtrack_idx].type = 0;  // CODE_END equivalent

    return backtrack_idx + 1;
}

/**
 * WFA2-based alignment - drop-in replacement for DPD()
 *
 * @param s1 Reference sequence (int8_t encoded)
 * @param s2 Read sequence (int8_t encoded)
 * @param backtrack Output backtrack array
 * @param size_neighbour_in_symbols Sequence length
 * @return Length of backtrack array, or -1 on error
 */
int DPD_WFA2(int8_t *s1, int8_t *s2, backtrack_t *backtrack, int size_neighbour_in_symbols) {
    // Convert to ASCII strings for WFA2
    char *pattern = (char *)malloc(size_neighbour_in_symbols + 1);
    char *text = (char *)malloc(size_neighbour_in_symbols + 1);

    if (!pattern || !text) {
        free(pattern);
        free(text);
        return -1;
    }

    upvc_to_ascii(s1, size_neighbour_in_symbols, pattern);
    upvc_to_ascii(s2, size_neighbour_in_symbols, text);

    // Configure WFA2 with UPVC's alignment parameters
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine;
    attributes.affine_penalties.match = 0;
    attributes.affine_penalties.mismatch = COST_SUB;
    attributes.affine_penalties.gap_opening = COST_GAPO;
    attributes.affine_penalties.gap_extension = COST_GAPE;

    // Use ultra-low memory mode (BiWFA) - O(s) space instead of O(n²)
    attributes.memory_mode = wavefront_memory_ultralow;

    // Enable alignment traceback
    attributes.alignment_scope = compute_alignment;

    // Initialize aligner
    wavefront_aligner_t *wf_aligner = wavefront_aligner_new(&attributes);
    if (!wf_aligner) {
        free(pattern);
        free(text);
        return -1;
    }

    // Perform alignment
    int status = wavefront_align(wf_aligner, pattern, size_neighbour_in_symbols,
                                  text, size_neighbour_in_symbols);

    if (status != 0) {
        // Alignment failed
        wavefront_aligner_delete(wf_aligner);
        free(pattern);
        free(text);
        return -1;
    }

    // Convert CIGAR to backtrack format
    int backtrack_len = cigar_to_backtrack(wf_aligner->cigar, s1, s2,
                                            backtrack, size_neighbour_in_symbols);

    // Cleanup
    wavefront_aligner_delete(wf_aligner);
    free(pattern);
    free(text);

    return backtrack_len;
}

/**
 * Optimized version that reuses aligner across multiple calls
 * Avoids allocation overhead for repeated alignments
 */
struct wfa2_aligner_pool_t {
    wavefront_aligner_t *aligner;
    char *pattern_buffer;
    char *text_buffer;
    int buffer_size;
};

wfa2_aligner_pool_t* wfa2_aligner_pool_new(int max_seq_len) {
    wfa2_aligner_pool_t *pool = (wfa2_aligner_pool_t *)malloc(sizeof(wfa2_aligner_pool_t));
    if (!pool) return NULL;

    // Configure aligner
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine;
    attributes.affine_penalties.match = 0;
    attributes.affine_penalties.mismatch = COST_SUB;
    attributes.affine_penalties.gap_opening = COST_GAPO;
    attributes.affine_penalties.gap_extension = COST_GAPE;
    attributes.memory_mode = wavefront_memory_ultralow;
    attributes.alignment_scope = compute_alignment;

    pool->aligner = wavefront_aligner_new(&attributes);
    pool->pattern_buffer = (char *)malloc(max_seq_len + 1);
    pool->text_buffer = (char *)malloc(max_seq_len + 1);
    pool->buffer_size = max_seq_len;

    if (!pool->aligner || !pool->pattern_buffer || !pool->text_buffer) {
        wfa2_aligner_pool_delete(pool);
        return NULL;
    }

    return pool;
}

void wfa2_aligner_pool_delete(wfa2_aligner_pool_t *pool) {
    if (!pool) return;
    if (pool->aligner) wavefront_aligner_delete(pool->aligner);
    free(pool->pattern_buffer);
    free(pool->text_buffer);
    free(pool);
}

int wfa2_aligner_pool_align(wfa2_aligner_pool_t *pool, int8_t *s1, int8_t *s2,
                             backtrack_t *backtrack, int seq_len) {
    if (!pool || seq_len > pool->buffer_size) return -1;

    // Convert sequences
    upvc_to_ascii(s1, seq_len, pool->pattern_buffer);
    upvc_to_ascii(s2, seq_len, pool->text_buffer);

    // Align
    int status = wavefront_align(pool->aligner, pool->pattern_buffer, seq_len,
                                  pool->text_buffer, seq_len);

    if (status != 0) return -1;

    // Convert CIGAR
    return cigar_to_backtrack(pool->aligner->cigar, s1, s2, backtrack, seq_len);
}
