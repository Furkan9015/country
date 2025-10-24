/**
 * WFA2 Wavefront Alignment Wrapper for UPVC
 * Header file
 */

#ifndef __WFA2_WRAPPER_H__
#define __WFA2_WRAPPER_H__

#include <stdint.h>
#include "upvc.h"  // For backtrack_t definition

/**
 * WFA2-based alignment - drop-in replacement for DPD()
 *
 * @param s1 Reference sequence (2-bit encoded)
 * @param s2 Read sequence (2-bit encoded)
 * @param backtrack Output backtrack array
 * @param size_neighbour_in_symbols Sequence length
 * @return Length of backtrack array, or -1 on error
 */
int DPD_WFA2(int8_t *s1, int8_t *s2, backtrack_t *backtrack, int size_neighbour_in_symbols);

/**
 * Reusable aligner pool for better performance
 */
typedef struct wfa2_aligner_pool_t wfa2_aligner_pool_t;

/**
 * Create a reusable aligner pool
 *
 * @param max_seq_len Maximum sequence length
 * @return Pool handle, or NULL on error
 */
wfa2_aligner_pool_t* wfa2_aligner_pool_new(int max_seq_len);

/**
 * Delete aligner pool
 */
void wfa2_aligner_pool_delete(wfa2_aligner_pool_t *pool);

/**
 * Align using pool (faster for repeated calls)
 *
 * @param pool Aligner pool
 * @param s1 Reference sequence
 * @param s2 Read sequence
 * @param backtrack Output backtrack array
 * @param seq_len Sequence length
 * @return Length of backtrack array, or -1 on error
 */
int wfa2_aligner_pool_align(wfa2_aligner_pool_t *pool, int8_t *s1, int8_t *s2,
                             backtrack_t *backtrack, int seq_len);

#endif /* __WFA2_WRAPPER_H__ */
