/**
 * OCOCO Counter-Based Variant Calling for UPVC
 *
 * This module replaces the traditional variant tree with OCOCO's
 * compact counter-based statistics for online/streaming variant calling.
 *
 * Stage 1: Basic Integration - Speed Optimization
 */

#ifndef __OCOCO_COUNTERS_H__
#define __OCOCO_COUNTERS_H__

#include <stdint.h>
#include <stdbool.h>
#include <pthread.h>

/**
 * Counter configuration:
 * - 16-bit mode: 4 bits per counter (max count 15), 4 counters = 16 bits total
 * - 32-bit mode: 7 bits per counter (max count 127), 4 counters + 4-bit ref = 32 bits
 * - 64-bit mode: 15 bits per counter (max count 32767), 4 counters + 4-bit ref = 64 bits
 *
 * We use 32-bit mode as a balance between memory and precision.
 */

#define COUNTER_BITS 7
#define COUNTER_MAX ((1 << COUNTER_BITS) - 1)
#define REFBASE_BITS 4

/**
 * Position statistics - compressed format (32 bits total)
 * Layout: [7 bits A][7 bits C][7 bits G][7 bits T][4 bits ref_nt16]
 */
typedef uint32_t pos_stat_compressed_t;

/**
 * Position statistics - uncompressed format for computation
 */
typedef struct {
    uint16_t counters[4];  // A, C, G, T counters
    uint8_t ref_nt16;      // Reference nucleotide (nt16 encoding)
    uint32_t sum;          // Total count (for quick access)
} pos_stat_uncompressed_t;

/**
 * OPTIMIZATION: Direct-access uncompressed counters (4 bytes per position)
 * Used during processing to eliminate 553M compress/decompress cycles
 *
 * Memory: 4 bytes × 78M positions × 8 threads = 2.5 GB (SAME as compressed!)
 * Speedup: Eliminates ~20 cycles per operation (decompress + compress overhead)
 * Max coverage: 255 per base (sufficient for typical 30-100x coverage)
 *
 * Trade-off: ZERO extra memory, potential speedup from direct access + better cache
 */
typedef struct {
    uint8_t counters[4];  // A, C, T, G (8-bit each, max 255)
} pos_stat_direct_t;

/**
 * Nucleotide encoding (nt4)
 * Maps A, C, T, G to 0, 1, 2, 3
 */
typedef enum {
    NT4_A = 0,
    NT4_C = 1,
    NT4_T = 2,
    NT4_G = 3,
    NT4_N = 4   // Unknown/ambiguous
} nt4_t;

/**
 * Variant entry for hash table (fixed-size, no malloc overhead)
 * Most variants are SNPs (1 char ref/alt), small indels fit in 15 chars
 */
#define MAX_ALLELE_LEN 15
typedef struct ococo_variant_entry {
    char ref[MAX_ALLELE_LEN + 1];           // Reference allele (inline, no malloc)
    char alt[MAX_ALLELE_LEN + 1];           // Alternate allele (inline, no malloc)
    _Atomic uint32_t depth;                 // Atomic depth counter
    _Atomic uint32_t score;                 // Atomic score accumulator
    _Atomic uint8_t occupied;               // Slot occupied flag
} ococo_variant_entry_t;

/**
 * Variant hash table per position (4 slots, handles SNPs + small indels)
 * Linear probing for collisions, rarely needed for typical variant density
 */
#define VARIANTS_PER_POSITION 4
typedef struct ococo_variant_table {
    ococo_variant_entry_t entries[VARIANTS_PER_POSITION];
} ococo_variant_table_t;

/**
 * Position list entry for fast variant iteration (Stage 2.5+ optimization)
 * Instead of scanning all genomic positions, we track only positions with variants
 */
typedef struct variant_position {
    uint32_t seq_nr;                          // Sequence (chromosome) number
    uint64_t pos;                             // Position within sequence
} variant_position_t;

/**
 * Dynamic position list for tracking variant positions
 * Used in OCOCO_HYBRID_MODE to avoid scanning entire genome
 */
typedef struct {
    variant_position_t *positions;            // Dynamic array of positions
    _Atomic uint64_t count;                   // Number of positions tracked
    uint64_t capacity;                        // Allocated capacity
    pthread_mutex_t lock;                     // Mutex for thread-safe insertions
} variant_position_list_t;

/**
 * Sparse variant table hash map entry
 * Only allocates variant tables for positions that actually have variants
 * Key: (seq_nr << 32) | pos
 */
typedef struct sparse_variant_entry {
    uint64_t key;                             // (seq_nr << 32) | pos
    ococo_variant_table_t *table;             // Variant hash table for this position
    struct sparse_variant_entry *next;       // Chaining for collisions
} sparse_variant_entry_t;

#define SPARSE_HASH_SIZE 8388608              // 8M buckets (reduced load factor from 70 to 10)

/**
 * Sparse variant table hash map
 * Memory: ~100K variants × 280 bytes = ~28 MB (vs 21.8 GB dense)
 */
typedef struct {
    sparse_variant_entry_t **buckets;         // Array of hash buckets
    pthread_mutex_t *locks;                   // Per-bucket locks for thread safety
    _Atomic uint64_t entry_count;             // Total number of entries
} sparse_variant_map_t;

/**
 * ATOMIC SPARSE COUNTER OPTIMIZATION (addresses memory-bound bottleneck)
 *
 * Three optimizations combined:
 * 1. ATOMIC: Single copy with atomic ops (vs 8 per-thread copies)
 * 2. SPARSE: Only allocate at positions with coverage (vs all 78M positions)
 * 3. 16-BIT: Direct 16-bit counters (vs 7-bit compressed), sufficient for 65K reads/position
 *
 * Memory savings:
 * - Old: 8 threads × 78M positions × 4 bytes = 2.5 GB
 * - New: ~45M covered positions × 9 bytes = ~405 MB
 * - Savings: 2.1 GB (84% reduction!)
 *
 * Performance improvement:
 * - Less memory bandwidth (single copy)
 * - Smaller working set (better cache hit rate)
 * - Atomic ops faster than memory bandwidth bottleneck
 */

/**
 * Atomic position counters (16-bit, uncompressed)
 * Sufficient for 30x average coverage (max 65535 per position)
 */
typedef struct {
    _Atomic uint16_t counters[4];  // A, C, T, G (atomic 16-bit each)
    uint8_t ref_base;               // Reference nucleotide (nt4 encoding)
} atomic_pos_counters_t;

/**
 * Sparse counter hash map entry
 * Only allocated for positions with actual coverage
 */
typedef struct sparse_counter_entry {
    uint64_t key;                           // (seq_nr << 32) | pos
    atomic_pos_counters_t *counters;        // Atomic counter array
    struct sparse_counter_entry *next;      // Chaining for collisions
} sparse_counter_entry_t;

#define SPARSE_COUNTER_HASH_SIZE 4194304    // 4M buckets for ~45M entries (load factor ~11)

/**
 * Sparse atomic counter map
 * Thread-safe atomic operations, on-demand allocation
 */
typedef struct {
    sparse_counter_entry_t **buckets;       // Array of hash buckets
    pthread_mutex_t *locks;                 // Per-bucket locks for get-or-create
    _Atomic uint64_t entry_count;           // Total positions with coverage
} sparse_counter_map_t;

/**
 * ========================================================================
 * DENSE ARRAY COUNTERS - ELIMINATES HASH TABLE SATURATION
 * ========================================================================
 *
 * CRITICAL FIX for DPU_OFFSET 224 slowdown:
 * - Hash tables accumulate entries across ALL batches (0-224)
 * - By batch 224: 87% full, average chain length 2.5, some chains 20-40 entries
 * - Each counter access requires O(n) chain traversal = exponential slowdown
 *
 * Dense array solution:
 * - Pre-allocate full genome-sized arrays (chr18: 82M positions)
 * - Direct array indexing: O(1) access, NO degradation over time
 * - Memory: 82M × 4 bytes × 8 threads = 2.6 GB (SAME as sparse with 58x coverage!)
 * - Eliminates: hash computation, collision resolution, malloc() calls, pointer chasing
 *
 * With 58x coverage, we access ~70-80M of 82M positions anyway, so "sparse"
 * optimization saves NO memory but adds massive overhead!
 */

/**
 * Per-thread dense counter array (NO hash tables, NO malloc during processing!)
 * Each processing thread has one of these (no synchronization needed!)
 */
typedef struct {
    pos_stat_direct_t *counters;            // Dense array: genome_size positions
    uint64_t size;                          // Number of positions (genome length)
} perthread_dense_array_t;

/**
 * Global OCOCO data structures
 *
 * DENSE ARRAY APPROACH (Fixes DPU_OFFSET 224 hash saturation):
 * - counters: PER-THREAD DENSE ARRAYS (full genome size, O(1) access)
 * - variants: Per-variant tracking (sparse hash map - only ~100K variants)
 *
 * Benefits:
 * - Speed: Constant O(1) access, NO progressive degradation
 * - Memory: Same as sparse (2.6 GB) since 58x coverage touches most positions
 * - Simplicity: Direct array indexing, no hash computation or collision handling
 *
 * Counter access: Direct array index = genome_position
 */
typedef struct {
    uint32_t nb_seq;                          // Number of sequences (chromosomes)
    uint64_t *len_seq;                        // Length of each sequence
    uint32_t num_threads;                     // Number of processing threads

    // PER-THREAD DENSE ARRAYS - pre-allocated for full genome
    // Each thread has its own array (no contention, O(1) access)
    perthread_dense_array_t *dense_counters;  // Array of per-thread dense arrays
    uint8_t **ref_bases;                      // Reference base at each position

    // Per-variant tracking (for accuracy) - SPARSE HASH MAP
    sparse_variant_map_t *variant_map;        // Sparse variant hash map (on-demand allocation)

    // STAGE 2.5+: Position list for fast variant iteration (OCOCO_HYBRID_MODE only)
    variant_position_list_t *position_list;   // Tracks positions with variants

    // Configuration parameters
    float majority_threshold;                 // Threshold for calling variant (default 0.5)
    uint32_t min_coverage;                    // Minimum coverage to call variant
    uint32_t init_ref_weight;                 // Initial weight for reference base
} ococo_counters_t;

/**
 * Initialize OCOCO counter arrays for the genome
 * Allocates memory for all genomic positions and all threads
 * @param num_threads Number of processing threads (for per-thread counters)
 */
void ococo_counters_init(uint32_t num_threads);

/**
 * Initialize reference bases from genome (call AFTER genome_load)
 */
void ococo_init_reference_bases(void);

/**
 * Free OCOCO counter arrays
 */
void ococo_counters_free(void);

/**
 * Get the global OCOCO counters structure
 */
ococo_counters_t *ococo_counters_get(void);

/**
 * Convert character nucleotide to nt4 encoding
 * @param nuc Character nucleotide (A, C, G, T, or other)
 * @return nt4 encoding (0-4)
 */
nt4_t char_to_nt4(char nuc);

/**
 * Convert nt4 encoding to character
 * @param nt4 Nucleotide encoding (0-4)
 * @return Character representation
 */
char nt4_to_char(nt4_t nt4);

/**
 * Increment counter for a specific nucleotide at a genomic position
 * DENSE ARRAY VERSION:
 * - Direct array access: counters[offset].counters[nt4]++
 * - O(1) constant time, no hash computation or collision handling
 * - No malloc during processing, perfect cache locality
 *
 * @param thread_id Thread ID for per-thread dense array access
 * @param seq_nr Sequence (chromosome) number
 * @param offset Position offset within the sequence
 * @param nt4 Nucleotide to increment (0-3 for A,C,T,G)
 */
void ococo_increment_counter(uint32_t thread_id, uint32_t seq_nr, uint64_t offset, nt4_t nt4);

/**
 * Decompress position statistics for analysis
 * @param compressed Compressed 32-bit counter data
 * @param uncompressed Output uncompressed statistics
 */
void ococo_decompress_stats(pos_stat_compressed_t compressed, pos_stat_uncompressed_t *uncompressed);

/**
 * Compress position statistics for storage
 * @param uncompressed Uncompressed statistics
 * @return Compressed 32-bit counter data
 */
pos_stat_compressed_t ococo_compress_stats(const pos_stat_uncompressed_t *uncompressed);

/**
 * Call consensus for a single position using majority rule
 * @param stats Uncompressed position statistics
 * @param threshold Majority threshold (typically 0.5)
 * @return Called nucleotide character
 */
char ococo_call_consensus_position(const pos_stat_uncompressed_t *stats, float threshold);

/**
 * Generate VCF file from OCOCO counters
 * This replaces the variant tree-based VCF generation
 *
 * @param vcf_filename Output VCF file path
 */
void ococo_create_vcf(const char *vcf_filename);

/**
 * Insert or update a variant in the lock-free variant list
 * Uses compare-and-swap for thread-safe insertion without mutex
 *
 * @param seq_nr Sequence (chromosome) number
 * @param pos Position within sequence
 * @param ref Reference allele string
 * @param alt Alternate allele string
 * @param score Smith-Waterman score for this variant instance
 */
void ococo_insert_variant(
    uint32_t thread_id,
    uint32_t seq_nr,
    uint64_t pos,
    const char *ref,
    const char *alt,
    uint32_t score
);

/**
 * Process alignment codes to extract and track variants
 * Parses CODE_SUB, CODE_INS, CODE_DEL from Smith-Waterman alignment
 *
 * HYBRID APPROACH (Per-Thread):
 * - Updates per-thread per-base counters (NO CAS overhead!)
 * - Extracts and inserts variants (for accuracy)
 *
 * @param thread_id Thread ID for per-thread counter access
 * @param seq_nr Sequence (chromosome) number
 * @param genome_pos Starting position in the genome
 * @param read Aligned read sequence
 * @param ref Reference genome sequence at alignment position
 * @param alignment_code Alignment code with SUB/INS/DEL operations
 * @param code_len Length of alignment code
 * @param size_read Length of the read
 */
void ococo_process_alignment(
    uint32_t thread_id,
    uint32_t seq_nr,
    uint64_t genome_pos,
    int8_t *read,
    int8_t *ref,
    uint8_t *alignment_code,
    int code_len,
    int size_read
);

/**
 * Print statistics about OCOCO counters
 * For debugging and performance analysis
 */
void ococo_print_stats(void);

#ifdef OCOCO_BAYESIAN_CALLING
/**
 * ========================================================================
 * BAYESIAN GENOTYPE CALLING (Phase 2+3)
 * ========================================================================
 */

/**
 * Genotype types for diploid calling
 */
typedef enum {
    GT_HOM_REF,   // 0/0 - Homozygous reference
    GT_HET,       // 0/1 - Heterozygous
    GT_HOM_ALT    // 1/1 - Homozygous alternate
} genotype_t;

/**
 * Genotype likelihoods (posterior probabilities)
 */
typedef struct {
    double hom_ref;   // P(0/0 | data)
    double het;       // P(0/1 | data)
    double hom_alt;   // P(1/1 | data)
} genotype_likelihoods_t;

/**
 * Genotype call result
 */
typedef struct {
    genotype_t genotype;    // Called genotype (0/0, 0/1, or 1/1)
    double confidence;      // Posterior probability of called genotype
    int gq;                 // Genotype quality (Phred-scaled)
} genotype_call_t;

/**
 * Calculate binomial log-likelihood: log P(k | n, p)
 * Uses log-gamma for numerical stability
 *
 * @param k Number of successes (e.g., alt allele count)
 * @param n Number of trials (e.g., total depth)
 * @param p Probability of success (e.g., 0.5 for het)
 * @return Log-likelihood value
 */
double binom_log_likelihood(int k, int n, double p);

/**
 * Calculate genotype likelihoods using Bayesian model
 *
 * @param ref_count Reference allele count
 * @param alt_count Alternate allele count
 * @param error_rate Sequencing error rate (typically 0.01)
 * @return Genotype likelihoods structure
 */
genotype_likelihoods_t calculate_genotype_likelihoods(
    uint16_t ref_count,
    uint16_t alt_count,
    double error_rate
);

/**
 * Call genotype using Bayesian inference
 * Returns the most likely genotype with confidence scores
 *
 * @param ref_count Reference allele count
 * @param alt_count Alternate allele count
 * @param error_rate Sequencing error rate
 * @return Genotype call with quality score
 */
genotype_call_t call_genotype_bayesian(
    uint16_t ref_count,
    uint16_t alt_count,
    double error_rate
);

/**
 * Get coverage-stratified error rate
 * Adapts error rate based on coverage depth
 *
 * @param cov Coverage depth
 * @return Estimated error rate for this coverage
 */
double get_error_rate_for_coverage(uint32_t cov);

/**
 * Convert genotype to VCF string format
 *
 * @param gt Genotype enum
 * @return VCF genotype string (e.g., "0/1")
 */
const char *genotype_to_vcf_string(genotype_t gt);

#endif /* OCOCO_BAYESIAN_CALLING */

#endif /* __OCOCO_COUNTERS_H__ */
