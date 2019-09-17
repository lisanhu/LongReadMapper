#ifndef _CHAINING_
#define _CHAINING_

#ifdef __cplusplus
#include <cstdint>
#else
#include <stdint.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
/// location type
typedef uint64_t ch_loc_t;
#define CH_LOC_MAX INT64_MAX

/// chain score type
typedef int32_t ch_score_t;
#define CH_SCORE_MIN INT32_MIN

/// chain size type
typedef uint32_t ch_size_t;

/**
 *  Anchor data structure
 *  It's a triple of x, y, w
 *      where x is the start location of the seed in the reference genome
 *          y is the start location of the seed in the read
 *          w is the seed length
 */
typedef struct ch_anchor_t {
    ch_loc_t x, y, w;
} ch_anchor_t;

/**
 *  A chain of anchors
 *  The score is for the chain
 */
typedef struct ch_chain_t {
    ch_score_t score;
    ch_size_t len;
    ch_loc_t pos;
//    ch_anchor_t anchors[];
} ch_chain_t;

/**
 *  The result of chaining
 *  The scores here are updated scores based on secondary chain
 */
typedef struct ch_result_t {
    ch_score_t scores[2];
    ch_chain_t chains[2];
} ch_result_t;

/**
 *  Chaining these anchors
 *  Each seed will produce several anchors
 *      with same y, and w
 *  The seed afterward should follow the previous seeds
 *      colinear
 *      non-overlapping
 */
ch_result_t chaining_seed_anchors(ch_anchor_t *_anchors,
                                  ch_loc_t num_anchros_per_seed,
                                  ch_loc_t num_seeds);

/**
 *  Chaining the anchors
 */
ch_result_t chaining(ch_anchor_t *anchors, ch_size_t num_anchors);

#ifdef __cplusplus
}
#endif
#endif
