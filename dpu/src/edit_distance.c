/**
 * Landau-Vishkin edit distance algorithm for DPU
 *
 * Computes edit distance in O(k*n) time using wavefront approach
 * Memory: O(k) instead of O(n^2) - perfect for DPU constraints
 *
 * Reference: "Fast text searching: allowing errors" (Landau & Vishkin, 1989)
 */

#include <stdint.h>
#include <string.h>
#include "edit_distance.h"

/* Maximum edit distance we care about - matches DPU MAX_SCORE threshold */
#define MAX_EDIT_DIST 40

/**
 * Landau-Vishkin edit distance algorithm
 *
 * @param s1 First sequence (reference neighbor)
 * @param s2 Second sequence (read neighbor)
 * @param max_dist Maximum distance to compute (early exit if exceeded)
 * @param len Length of sequences in bytes
 * @return Edit distance, or UINT32_MAX if exceeds max_dist
 */
uint32_t edit_distance_landau_vishkin(uint8_t *s1, uint8_t *s2, uint32_t max_dist, uint32_t len)
{
    /* Wavefront array: for each diagonal d at distance k, stores furthest row reached
     * We need 2*max_dist+1 diagonals centered around main diagonal */
    int32_t farthest[2 * MAX_EDIT_DIST + 1];
    int32_t prev_farthest[2 * MAX_EDIT_DIST + 1];

    /* Limit max distance to our buffer size */
    if (max_dist > MAX_EDIT_DIST) {
        max_dist = MAX_EDIT_DIST;
    }

    /* Center index for diagonal 0 */
    const int32_t center = (int32_t)max_dist;

    /* Initialize: diagonal 0 at distance 0 starts at position 0 */
    for (int32_t d = -center; d <= center; d++) {
        prev_farthest[d + center] = -1;
    }
    prev_farthest[center] = 0;

    /* Match as far as possible from position 0 on main diagonal */
    int32_t row = 0;
    while (row < (int32_t)len && (s1[row] & 0x3) == (s2[row] & 0x3)) {
        row++;
    }
    prev_farthest[center] = row;

    /* Check if perfect match */
    if (row == (int32_t)len) {
        return 0;
    }

    /* Iterate over edit distances k = 1, 2, 3, ... */
    for (uint32_t k = 1; k <= max_dist; k++) {
        /* For each diagonal d in range [-k, +k] */
        for (int32_t d = -(int32_t)k; d <= (int32_t)k; d++) {
            const int32_t diag_idx = d + center;
            int32_t max_row = -1;

            /* Three ways to reach diagonal d at distance k:
             * 1. Substitution: from diagonal d at distance k-1
             * 2. Deletion: from diagonal d-1 at distance k-1, then move right
             * 3. Insertion: from diagonal d+1 at distance k-1, then move down
             */

            /* Substitution/Match: same diagonal, previous distance */
            if (d >= -(int32_t)(k-1) && d <= (int32_t)(k-1)) {
                int32_t r = prev_farthest[diag_idx] + 1;  /* Move to next row */
                if (r > max_row) max_row = r;
            }

            /* Deletion: from left diagonal */
            if (d > -(int32_t)k && diag_idx - 1 >= 0) {
                int32_t r = prev_farthest[diag_idx - 1] + 1;  /* Move right */
                if (r > max_row) max_row = r;
            }

            /* Insertion: from right diagonal */
            if (d < (int32_t)k && diag_idx + 1 <= 2 * center) {
                int32_t r = prev_farthest[diag_idx + 1];  /* Move down */
                if (r > max_row) max_row = r;
            }

            /* Extend matches along this diagonal as far as possible */
            int32_t row = max_row;
            int32_t col = max_row + d;
            while (row < (int32_t)len && col < (int32_t)len && col >= 0 &&
                   (s1[row] & 0x3) == (s2[col] & 0x3)) {
                row++;
                col++;
            }

            farthest[diag_idx] = row;

            /* Check if we reached the end */
            if (row >= (int32_t)len && col >= (int32_t)len) {
                return k;
            }
        }

        /* Swap buffers for next iteration */
        for (int32_t d = -(int32_t)k; d <= (int32_t)k; d++) {
            prev_farthest[d + center] = farthest[d + center];
        }
    }

    /* Distance exceeds max_dist */
    return UINT32_MAX;
}
