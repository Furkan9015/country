/**
 * Edit distance computation for DPU
 */

#ifndef __EDIT_DISTANCE_H__
#define __EDIT_DISTANCE_H__

#include <stdint.h>

/**
 * @brief Computes edit distance using Landau-Vishkin algorithm
 *
 * @param s1 First sequence (reference neighbor)
 * @param s2 Second sequence (read neighbor)
 * @param max_dist Maximum edit distance to compute
 * @param len Length of sequences in bytes
 * @return Edit distance, or UINT32_MAX if exceeds max_dist
 */
uint32_t edit_distance_landau_vishkin(uint8_t *s1, uint8_t *s2, uint32_t max_dist, uint32_t len);

#endif /* __EDIT_DISTANCE_H__ */
