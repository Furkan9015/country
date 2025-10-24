/**
 * @file mram_counters.h
 * @brief DPU MRAM-based variant counters for PIM acceleration
 *
 * Stage 2 of OCOCO-PIM architecture: Move counter updates to DPU MRAM
 * to eliminate host processing bottleneck.
 */

#ifndef MRAM_COUNTERS_H
#define MRAM_COUNTERS_H

#include <stdint.h>
#include <mram.h>

/**
 * @brief Per-position counter structure in MRAM
 *
 * Stores counts for each nucleotide base (A, C, G, T) at a genomic position.
 * Total size: 8 bytes (aligned for DMA - DMA transfers must be 8-byte aligned)
 */
typedef struct {
    uint8_t A;  ///< Counter for adenine
    uint8_t C;  ///< Counter for cytosine
    uint8_t G;  ///< Counter for guanine
    uint8_t T;  ///< Counter for thymine
    uint8_t reserved[4];  ///< Padding for 8-byte DMA alignment
} __attribute__((aligned(8))) mram_counter_t;

/**
 * @brief Initialize DPU counter array in MRAM
 *
 * Allocates and zero-initializes counter array for this DPU's genome slice.
 * Called once during DPU initialization.
 *
 * @param slice_length Number of bases in this DPU's genome slice
 */
void mram_counters_init(uint32_t slice_length);

/**
 * @brief Update counters based on alignment CIGAR
 *
 * Processes CIGAR operations to increment base counters at aligned positions.
 * Core function for Stage 2 - moves variant counting to DPU.
 *
 * @param cigar_ops     Array of CIGAR operations (MATCH/INS/DEL)
 * @param cigar_len     Number of CIGAR operations
 * @param ref_offset    Starting position in reference slice
 * @param read_seq      Read sequence (2-bit packed)
 * @param stats         Statistics structure to update
 */
void mram_update_counters_from_cigar(
    const uint8_t *cigar_ops,
    uint8_t cigar_len,
    uint32_t ref_offset,
    const uint8_t *read_seq,
    void *stats
);

/**
 * @brief Get base pointer to counter array in MRAM
 *
 * @return Pointer to start of counter array in MRAM
 */
__mram_ptr mram_counter_t *mram_get_counter_base(void);

/**
 * @brief Get size of counter array in bytes
 *
 * @return Size of counter array (slice_length * sizeof(mram_counter_t))
 */
uint32_t mram_get_counter_size(void);

#endif /* MRAM_COUNTERS_H */
