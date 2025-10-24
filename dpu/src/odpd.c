/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include "odpd.h"
#include "debug.h"
#include "common.h"  /* For CIGAR_OP_* constants */
#include <alloc.h>
#include <defs.h>    /* For me() function */
#include <mram.h>    /* For MRAM backtrack storage */

#define NB_DIAGS 15
#define COST_SUB 10
#define COST_GAPO 11
#define COST_GAPE 1

static inline int min(int a, int b)
{
    if (a < b)
        return a;
    else
        return b;
}

/* Need three huge tables per tasklet to work. */
static int *__D;
static int *__P;
static int *__Q;

/* MRAM storage for backtrack matrices (1 byte per cell, shared among all tasklets) */
__mram_noinit uint8_t m_backtrack_storage[NR_TASKLETS][110][110];

/* Macros for MRAM backtrack access with WRAM caching */
#define BACKTRACK_ROW_SIZE 112  /* 110 rounded up to 8-byte boundary for DMA */

static inline unsigned int nb_items_per_matrix(unsigned int nbr_sym_len)
{
    /* One matrix is ...[2][nbr_sym_len + 1] */
    return (nbr_sym_len + 1) << 1;
}

static inline unsigned int sizeof_matrix(unsigned int nbr_sym_len)
{
    /* Each matrix contains ints = 4 bytes */
    return nb_items_per_matrix(nbr_sym_len) << 2;
}

void odpd_init(unsigned int nb_tasklets, unsigned int nbr_sym_len)
{
    __D = mem_alloc(sizeof_matrix(nbr_sym_len) * nb_tasklets);
    __P = mem_alloc(sizeof_matrix(nbr_sym_len) * nb_tasklets);
    __Q = mem_alloc(sizeof_matrix(nbr_sym_len) * nb_tasklets);
}

static inline int *get_matrix_for_tasklet(int *base, unsigned int tid, unsigned int nbr_sym_len)
{
    uint8_t *as_bytes = (uint8_t *)base;
    uint8_t *result = as_bytes + sizeof_matrix(nbr_sym_len) * tid;
    return (int *)result;
}

static inline int *_at(int *M, int x, int y, int nbr_sym_len)
{
    /* x is equal to 0 or 1 */
    return M + (((x == 0) ? 0 : nbr_sym_len + 1) + y);
}

#define D(x, y) *_at(_D, x, y, nbr_sym_len)
#define P(x, y) *_at(_P, x, y, nbr_sym_len)
#define Q(x, y) *_at(_Q, x, y, nbr_sym_len)

/* C implementation of odpd removed - using optimized assembly version from odpd_opt.S */

/**
 * @brief Enhanced odpd with CIGAR output (Option B optimization)
 *
 * This version performs the same banded Smith-Waterman alignment as odpd(),
 * but also tracks the alignment path and outputs a CIGAR string.
 * Memory overhead: ~nbr_sym_len^2 bytes for backtrack matrix (allocated in WRAM per-tasklet).
 */
int odpd_with_cigar(const uint8_t *s1, const uint8_t *s2, int max_score, unsigned int nbr_sym_len,
                    uint8_t *cigar_out, uint8_t *cigar_len)
{
    unsigned int tid = me();
    int *_D = get_matrix_for_tasklet(__D, tid, nbr_sym_len);
    int *_P = get_matrix_for_tasklet(__P, tid, nbr_sym_len);
    int *_Q = get_matrix_for_tasklet(__Q, tid, nbr_sym_len);

    /* MRAM backtrack storage with WRAM row cache */
    __mram_ptr uint8_t *backtrack_mram = (__mram_ptr uint8_t *)&m_backtrack_storage[tid][0][0];
    __dma_aligned uint8_t backtrack_row_cache[BACKTRACK_ROW_SIZE];  /* Cache for current row */

    int i, j, d, lp, pp, QP, min_score;
    int final_i = nbr_sym_len, final_j = nbr_sym_len;

    /* Initialize DP matrices (same as odpd) */
    for (j = 0; j <= NB_DIAGS / 2 + 1; j++) {
        P(0, j) = 99;
        Q(0, j) = 99;
    }
    P(1, 0) = 99;
    Q(1, 0) = 99;
    for (j = 0; j <= NB_DIAGS / 2 + 1; j++) {
        D(0, j) = j * COST_SUB;
    }

    /* Forward DP with backtrack recording */
    for (i = 1; i < NB_DIAGS / 2 + 1; i++) {
        min_score = 99;
        pp = i & 1;
        lp = pp ^ 1;
        D(pp, 0) = i * COST_SUB;

        /* Clear row cache */
        for (int k = 0; k < BACKTRACK_ROW_SIZE; k++) {
            backtrack_row_cache[k] = 0;
        }

        for (j = 1; j <= i + NB_DIAGS / 2; j++) {
            P(pp, j) = min(D(pp, j - 1) + COST_GAPO, P(pp, j - 1) + COST_GAPE);
            Q(pp, j) = min(D(lp, j) + COST_GAPO, Q(lp, j) + COST_GAPE);
            QP = min(P(pp, j), Q(pp, j));
            d = D(lp, j - 1);
            int is_mismatch = (((s1[(i - 1) / 4] >> (2 * ((i - 1) % 4))) & 3) !=
                               ((s2[(j - 1) / 4] >> (2 * ((j - 1) % 4))) & 3));
            if (is_mismatch) {
                d += COST_SUB;
            }
            D(pp, j) = min(d, QP);

            /* Record backtrack decision in row cache */
            if (D(pp, j) == d) {
                backtrack_row_cache[j] = CIGAR_OP_MATCH;  /* Match or substitution */
            } else if (D(pp, j) == P(pp, j)) {
                backtrack_row_cache[j] = CIGAR_OP_INS;    /* Insertion */
            } else {
                backtrack_row_cache[j] = CIGAR_OP_DEL;    /* Deletion */
            }

            if (D(pp, j) < min_score) {
                min_score = D(pp, j);
            }
        }

        /* Write completed row to MRAM */
        mram_write(backtrack_row_cache, backtrack_mram + (i * 110), BACKTRACK_ROW_SIZE);

        Q(pp, j) = 99;
        D(pp, j) = 99;
        if (min_score > max_score) {
            *cigar_len = 0;  /* Score too high, don't return CIGAR */
            return min_score;
        }
    }

    /* Middle section */
    for (i = NB_DIAGS / 2 + 1; i < nbr_sym_len - NB_DIAGS / 2; i++) {
        min_score = 99;
        pp = i & 1;
        lp = pp ^ 1;
        j = i - NB_DIAGS / 2 - 1;
        P(pp, j) = 99;
        D(pp, j) = 99;

        /* Clear row cache */
        for (int k = 0; k < BACKTRACK_ROW_SIZE; k++) {
            backtrack_row_cache[k] = 0;
        }

        for (j = i - NB_DIAGS / 2; j <= i + NB_DIAGS / 2; j++) {
            P(pp, j) = min(D(pp, j - 1) + COST_GAPO, P(pp, j - 1) + COST_GAPE);
            Q(pp, j) = min(D(lp, j) + COST_GAPO, Q(lp, j) + COST_GAPE);
            QP = min(P(pp, j), Q(pp, j));
            d = D(lp, j - 1);
            int is_mismatch = (((s1[(i - 1) / 4] >> (2 * ((i - 1) % 4))) & 3) !=
                               ((s2[(j - 1) / 4] >> (2 * ((j - 1) % 4))) & 3));
            if (is_mismatch) {
                d += COST_SUB;
            }
            D(pp, j) = min(d, QP);

            /* Record backtrack in row cache */
            if (D(pp, j) == d) {
                backtrack_row_cache[j] = CIGAR_OP_MATCH;
            } else if (D(pp, j) == P(pp, j)) {
                backtrack_row_cache[j] = CIGAR_OP_INS;
            } else {
                backtrack_row_cache[j] = CIGAR_OP_DEL;
            }

            if (D(pp, j) < min_score) {
                min_score = D(pp, j);
            }
        }

        /* Write completed row to MRAM */
        mram_write(backtrack_row_cache, backtrack_mram + (i * 110), BACKTRACK_ROW_SIZE);

        Q(pp, j) = 99;
        D(pp, j) = 99;
        if (min_score > max_score) {
            *cigar_len = 0;
            return min_score;
        }
    }

    /* Final section */
    min_score = 99;
    for (i = nbr_sym_len - NB_DIAGS / 2; i < nbr_sym_len + 1; i++) {
        pp = i & 1;
        lp = pp ^ 1;
        j = i - NB_DIAGS / 2 - 1;
        P(pp, j) = 99;
        D(pp, j) = 99;

        /* Clear row cache */
        for (int k = 0; k < BACKTRACK_ROW_SIZE; k++) {
            backtrack_row_cache[k] = 0;
        }

        for (j = i - NB_DIAGS / 2; j < nbr_sym_len + 1; j++) {
            P(pp, j) = min(D(pp, j - 1) + COST_GAPO, P(pp, j - 1) + COST_GAPE);
            Q(pp, j) = min(D(lp, j) + COST_GAPO, Q(lp, j) + COST_GAPE);
            QP = min(P(pp, j), Q(pp, j));
            d = D(lp, j - 1);
            int is_mismatch = (((s1[(i - 1) / 4] >> (2 * ((i - 1) % 4))) & 3) !=
                               ((s2[(j - 1) / 4] >> (2 * ((j - 1) % 4))) & 3));
            if (is_mismatch) {
                d += COST_SUB;
            }
            D(pp, j) = min(d, QP);

            /* Record backtrack in row cache */
            if (D(pp, j) == d) {
                backtrack_row_cache[j] = CIGAR_OP_MATCH;
            } else if (D(pp, j) == P(pp, j)) {
                backtrack_row_cache[j] = CIGAR_OP_INS;
            } else {
                backtrack_row_cache[j] = CIGAR_OP_DEL;
            }
        }

        /* Write completed row to MRAM */
        mram_write(backtrack_row_cache, backtrack_mram + (i * 110), BACKTRACK_ROW_SIZE);

        if (D(pp, nbr_sym_len) < min_score) {
            min_score = D(pp, nbr_sym_len);
            final_i = i;
            final_j = nbr_sym_len;
        }
    }

    /* Find best ending position in last row */
    i = nbr_sym_len;
    pp = i & 1;
    for (j = i - NB_DIAGS / 2; j < nbr_sym_len + 1; j++) {
        if (D(pp, j) < min_score) {
            min_score = D(pp, j);
            final_i = i;
            final_j = j;
        }
    }

    /* Traceback from (final_i, final_j) to (0, 0) reading from MRAM */
    i = final_i;
    j = final_j;
    *cigar_len = 0;
    int cached_row = -1;  /* Track which row is currently cached */

    while (i > 0 && j > 0) {
        /* Load row from MRAM if not already cached */
        if (cached_row != i) {
            mram_read(backtrack_mram + (i * 110), backtrack_row_cache, BACKTRACK_ROW_SIZE);
            cached_row = i;
        }

        uint8_t op = backtrack_row_cache[j];
        cigar_out[*cigar_len] = op;
        (*cigar_len)++;

        if (op == CIGAR_OP_MATCH) {
            i--;
            j--;
        } else if (op == CIGAR_OP_INS) {
            j--;  /* Insertion in read */
        } else {  /* CIGAR_OP_DEL */
            i--;  /* Deletion in read */
        }

        if (*cigar_len >= MAX_CIGAR_OPS - 1) {
            break;  /* Safety: prevent overflow */
        }
    }

    /* CIGAR is in reverse order (from end to start), reverse it */
    for (uint8_t k = 0; k < *cigar_len / 2; k++) {
        uint8_t temp = cigar_out[k];
        cigar_out[k] = cigar_out[*cigar_len - 1 - k];
        cigar_out[*cigar_len - 1 - k] = temp;
    }

    return min_score;
}
