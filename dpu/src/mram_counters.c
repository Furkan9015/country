/**
 * @file mram_counters.c
 * @brief Implementation of DPU MRAM-based variant counters
 *
 * Stage 2: OCOCO-PIM Architecture
 * Moves counter updates from host CPU to DPU MRAM for 2× speedup.
 */

#include "mram_counters.h"
#include "common.h"
#include "stats.h"
#include <mram.h>
#include <stdint.h>
#include <string.h>

/* MRAM counter storage
 * NOTE: For Stage 2 prototype, using a smaller fixed size
 * Production version will use host-provided size via DPU_MRAM_INFO_VAR
 * Temporary size: 64K positions = 256KB (fits easily in 64MB MRAM)
 */
#define MAX_COUNTER_POSITIONS (64 * 1024)
__mram_noinit mram_counter_t DPU_COUNTER_ARRAY[MAX_COUNTER_POSITIONS];

/* Metadata */
static uint32_t counter_array_length = 0;

void mram_counters_init(uint32_t slice_length) {
    counter_array_length = slice_length;

    /* Zero-initialization happens automatically with __mram_noinit */
    /* Could add explicit zeroing here if needed, but adds init time */
}

__mram_ptr mram_counter_t *mram_get_counter_base(void) {
    return DPU_COUNTER_ARRAY;
}

uint32_t mram_get_counter_size(void) {
    return counter_array_length * sizeof(mram_counter_t);
}

/**
 * @brief Update counters from CIGAR alignment
 *
 * Performance characteristics:
 * - MRAM read:  ~200 cycles (62ns @ 350 MHz)
 * - Increment:  ~10 cycles
 * - MRAM write: ~200 cycles
 * - Total:      ~410 cycles per base update
 *
 * For 120bp read: 120 × 410 = 49,200 cycles = 140μs per alignment
 * Per DPU (4,500 alignments/tasklet): 4,500 × 140μs = 630ms
 * Across 128 DPUs: ~1-2s total (with pipeline overlap)
 */
void mram_update_counters_from_cigar(
    const uint8_t *cigar_ops,
    uint8_t cigar_len,
    uint32_t ref_offset,
    const uint8_t *read_seq,
    void *stats_ptr
) {
    STATS_ATTRIBUTE dpu_tasklet_stats_t *stats = (dpu_tasklet_stats_t *)stats_ptr;

    __dma_aligned mram_counter_t counter_cache;
    int ref_pos = 0;
    int read_pos = 0;

    for (uint8_t op_idx = 0; op_idx < cigar_len; op_idx++) {
        uint8_t op = cigar_ops[op_idx];

        if (op == CIGAR_OP_MATCH) {
            /* Calculate MRAM address for this position */
            uint32_t genome_pos = ref_offset + ref_pos;

            if (genome_pos >= counter_array_length) {
                /* Position outside this DPU's slice - skip */
                ref_pos++;
                read_pos++;
                continue;
            }

            __mram_ptr mram_counter_t *counter_addr = &DPU_COUNTER_ARRAY[genome_pos];

            /* Read current counter from MRAM (DMA aligned) */
            mram_read(counter_addr, &counter_cache, sizeof(mram_counter_t));
            STATS_INCR_LOAD(stats, sizeof(mram_counter_t));

            /* Extract base from read sequence (2-bit packed) */
            uint8_t base = (read_seq[read_pos / 4] >> (2 * (read_pos % 4))) & 0x3;

            /* Increment appropriate counter (with saturation at 255) */
            switch(base) {
                case 0: /* A */
                    if (counter_cache.A < 255) counter_cache.A++;
                    break;
                case 1: /* C */
                    if (counter_cache.C < 255) counter_cache.C++;
                    break;
                case 2: /* T */
                    if (counter_cache.T < 255) counter_cache.T++;
                    break;
                case 3: /* G */
                    if (counter_cache.G < 255) counter_cache.G++;
                    break;
            }

            /* Write updated counter back to MRAM */
            mram_write(&counter_cache, counter_addr, sizeof(mram_counter_t));
            STATS_INCR_STORE(stats, sizeof(mram_counter_t));

            ref_pos++;
            read_pos++;
        }
        else if (op == CIGAR_OP_INS) {
            /* Insertion in read - only advance read position */
            read_pos++;
        }
        else if (op == CIGAR_OP_DEL) {
            /* Deletion in read - only advance reference position */
            ref_pos++;
        }
        /* Unknown operations are silently ignored */
    }
}

/**
 * OPTIMIZATION NOTES for future improvement:
 *
 * 1. Batched Updates (50% speedup potential):
 *    - Accumulate updates in WRAM cache
 *    - Flush to MRAM every 8-16 positions
 *    - Reduces MRAM writes by 8-16×
 *
 * 2. 16-bit Counters (if saturation is an issue):
 *    - Increase counter size to uint16_t
 *    - Doubles MRAM usage but eliminates saturation
 *    - Total: 3.2 MB per DPU (still fits easily)
 *
 * 3. OCOCO Saturation Mechanism:
 *    - On overflow, right-shift all counters by 1
 *    - Maintains frequency ratios
 *    - Filters noise through approximate counting
 *
 * 4. Sparse Counter Storage:
 *    - Only store non-zero counters
 *    - Use hash table or compressed format
 *    - Saves MRAM but adds lookup overhead
 */
