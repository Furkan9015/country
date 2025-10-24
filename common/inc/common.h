/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __COMMON_H__
#define __COMMON_H__

#include <assert.h>
#include <stdint.h>

/* #define STATS_ON */

#define MRAM_SIZE (64 << 20)

#define ALIGN_DPU(val) (((val) + 7) & ~7)

#define MAX_DPU_REQUEST (1 << 15)
#define MAX_DPU_RESULTS (1 << 19)
#define MAX_RESULTS_PER_READ (1 << 10)

#define SIZE_READ 150
#define SIZE_SEED 14
#define SIZE_NEIGHBOUR_IN_BYTES ((SIZE_READ - SIZE_SEED) / 4)
#define DELTA_NEIGHBOUR(round) ((SIZE_SEED * round) / 4)
#define SIZE_IN_SYMBOLS(delta) ((SIZE_NEIGHBOUR_IN_BYTES - delta) * 4)

/**
 * @brief A snapshot of the MRAM.
 *
 * @var delta           Delta to apply to nodp and odpd comparison depending on the round.
 */
typedef uint32_t delta_info_t;
#define DPU_MRAM_INFO_VAR m_mram_info

/**
 * @brief Coordonates of the read that matched in the reference genome.
 */
typedef struct {
    union {
        uint64_t coord;
        struct {
            uint32_t seed_nr;
            uint32_t seq_nr;
        };
    };
} dpu_result_coord_t;

/**
 * @brief CIGAR operation encoding for DPU alignment results
 *
 * 2-bit encoding per operation:
 * 00 = Match/Mismatch (M)
 * 01 = Insertion (I)
 * 10 = Deletion (D)
 * 11 = Reserved
 */
#define CIGAR_OP_MATCH 0
#define CIGAR_OP_INS   1
#define CIGAR_OP_DEL   2

/**
 * @brief Maximum CIGAR operations for SIZE_READ (120bp)
 * Worst case: all mismatches/indels = 120 ops
 */
#define MAX_CIGAR_OPS 140

/**
 * @brief Bitpacked CIGAR storage
 * 2 bits per operation, 16 ops per uint32_t
 * 140 ops / 16 = 8.75, round up to 9 uint32_t = 36 bytes (vs 140 bytes unpacked)
 */
#define CIGAR_PACKED_WORDS 9

/**
 * @brief Helper macros for CIGAR bitpacking (2 bits per op)
 */
#define CIGAR_PACK_OP(packed, idx, op) \
    ((packed)[(idx) / 16] |= ((uint32_t)(op) << (((idx) % 16) * 2)))

#define CIGAR_UNPACK_OP(packed, idx) \
    (((packed)[(idx) / 16] >> (((idx) % 16) * 2)) & 0x3)

/**
 * @brief One result produced for one read
 *
 * @var num  Number that associate an input read with a request.
 * @var score Best score of the read with a read of the reference genome.
 * @var coord Coordinate of the read that matched in the reference genome.
 * @var cigar_len Number of CIGAR operations (0 = no CIGAR available)
 * @var cigar_packed Bitpacked CIGAR encoding (2 bits per operation, 16 ops per uint32_t)
 */
typedef struct {
    union {
        struct {
            uint32_t score;
            int32_t num;
        };
        uint64_t key;
    };
    dpu_result_coord_t coord;

    /* DPU CIGAR output (Option B optimization - bitpacked for memory efficiency) */
    uint8_t cigar_len;                       /* 0 = no CIGAR (fallback to host alignment) */
    uint8_t reserved[3];                     /* Reserved for future use */
    uint32_t cigar_packed[CIGAR_PACKED_WORDS]; /* Bitpacked CIGAR: 2 bits/op, 36 bytes */
} __attribute__((aligned(8))) dpu_result_out_t;
#define DPU_RESULT_VAR m_dpu_result

/**
 * @brief Number of results produces by all the tasklets
 */
typedef uint32_t nb_result_t;
#define DPU_NB_RESULT_VAR m_dpu_nb_result

/**
 * @brief stats reported by every tasklet
 */
typedef struct {
    uint32_t nb_reqs;
    uint32_t nb_nodp_calls;
    uint32_t nb_odpd_calls;
    uint32_t nb_results;
    uint32_t mram_data_load;
    uint32_t mram_result_store;
    uint32_t mram_load;
    uint32_t mram_store;
    uint64_t nodp_time;
    uint64_t odpd_time;
} dpu_tasklet_stats_t;
#define DPU_TASKLET_STATS_VAR m_dpu_tasklet_stats

typedef uint64_t dpu_compute_time_t;
#define DPU_COMPUTE_TIME_VAR m_dpu_compute_time

/**
 * @brief Information on the requests reads to a DPU.
 */
typedef uint32_t nb_request_t;
#define DPU_NB_REQUEST_VAR m_dpu_nb_request

/**
 * @brief Structure representing one request to a DPU, resulting from the dispatching of reads.
 *
 * Such a request basically contains all the information for one read.
 *
 * @var offset  The 1st neighbour address.
 * @var count   The number of neighbours.
 * @var num     A reference number to the original request.
 * @var nbr     The input neighbour to compare with the reference
 */
typedef struct {
    uint32_t offset;
    uint32_t count;
    uint32_t num;
    uint8_t nbr[SIZE_NEIGHBOUR_IN_BYTES];
} dpu_request_t;
#define DPU_REQUEST_VAR m_dpu_request

typedef struct {
    dpu_result_coord_t coord;
    uint8_t nbr[ALIGN_DPU(SIZE_NEIGHBOUR_IN_BYTES)];
} coords_and_nbr_t;

#endif /* __COMMON_H__ */
