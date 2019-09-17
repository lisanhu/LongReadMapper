#include "chain.h"
#include <cstdio>
#include <iostream>
#include "utils.h"

static ch_chain_t back_tracking(ch_size_t* ids, ch_size_t end);

ch_result_t chaining_seed_anchors(ch_anchor_t* _anchors,
                                  ch_loc_t num_anchros_per_seed,
                                  ch_loc_t num_seeds);
namespace chain {
ch_anchor_t const zero = {0, 0, 0};
ch_anchor_t const term = {CH_LOC_MAX, CH_LOC_MAX, 0};
ch_score_t const ge = -2;
ch_score_t const em = 2;

bool precede(ch_anchor_t const& a, ch_anchor_t const& b) {
    return a.y < b.y && a.x + a.w <= b.x;
}

bool operator==(ch_anchor_t const& lhs, ch_anchor_t const& rhs) {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.w == rhs.w;
}

ch_score_t gap(ch_anchor_t const& a, ch_anchor_t const& b) {
    if (!precede(a, b)) {
        return CH_SCORE_MIN;
    }
    if (a == zero || b == term) {
        return 0;
    }
    return ge * (b.x - a.x - a.w);
}

ch_score_t weight(ch_anchor_t const& a) { return em * a.w; }

ch_size_t find_precedes(ch_anchor_t a, ch_anchor_t const* anchors,
                        ch_size_t num_anchors, ch_anchor_t* buf) {
    ch_size_t len = 1;
    for (ch_size_t i = 0; i < num_anchors; ++i) {
        if (precede(anchors[i], a)) {
            buf[len++] = anchors[i];
        }
    }
    return len;
}

int anchor_to_string(ch_anchor_t const& anc, char* buf) {
    return sprintf(buf, "%lu Anchor{ x: %lu, y: %lu, w: %lu}", anc.x - anc.y,
                   anc.x, anc.y, anc.w);
}

}  // namespace chain

ch_anchor_t* ancs;
ch_chain_t back_tracking(ch_size_t* ids, ch_size_t end) {
    using namespace std;

    ch_size_t count = 0;
    printf("Start back_tracking...\n");
    char buf[100] = {};
    char buf2[100] = {};
    while (count++, ids[end - 1] != 0) {
        chain::anchor_to_string(ancs[end - 1], buf);
        chain::anchor_to_string(ancs[ids[end] - 1], buf2);
        printf("%3d: [%u, %u] %s ==> %s\n", count, end, ids[end - 1], buf, buf2);
        end = ids[end];
    }
    ch_chain_t chain = {0, count, end};
    return chain;
}

ch_result_t chaining(ch_anchor_t* anchors, ch_size_t num_anchors) {
    using namespace std;

    ancs = anchors;

    /// scores include zero and term
    ch_score_t scores[num_anchors + 2];
    ch_size_t ids[num_anchors + 2];

    for (ch_size_t i = 0; i < num_anchors; ++i) {
        scores[i + 1] = ids[i + 1] = 0;
    }
    /// score of zero is 0
    scores[0] = 0;
    for (ch_size_t i = 0; i < num_anchors; ++i) {
        ch_anchor_t& anc = anchors[i];
        // scores[i + 1] = chain::weight(anc);
        scores[i + 1] = 0;
        std::pair<int, std::size_t> max;
        typedef struct {
            ch_anchor_t* anchor;
        } Context;
        Context ctx;
        ctx.anchor = &anc;
        auto max_score = [](ch_anchor_t const& a, void const* context) -> int {
            Context const ctx = *static_cast<Context const*>(context);
            ch_score_t score = chain::weight(a) + chain::gap(a, *ctx.anchor);
            char buf[100];
            chain::anchor_to_string(*ctx.anchor, buf);
            printf("%s\n", buf);
            chain::anchor_to_string(a, buf);
            printf("%s\n", buf);
            // if (score > -2147483608) {
            printf("%lu, %lu => %d\n", a.y, ctx.anchor->y, score);
            //}
            return chain::weight(a) + chain::gap(a, *ctx.anchor);
        };
        // max = utils::arg_max(max_score, &ctx, anchors, num_anchors);
        max = utils::arg_max(max_score, &ctx, anchors, i);
        printf("id, max, max_id = %u, %d, %lu\n", i, max.first, max.second);
        printf("max, max_id = %d, %lu\n", max.first, max.second);
        if (max.first < 0) {
            max.first = 0, max.second = -1;
        }
        ids[i + 1] = max.second + 1;
        scores[i + 1] += max.first;
        scores[i + 1] += scores[max.second + 1];
        printf("score: %d, ", scores[i + 1]);
        printf("===========================\n");
    }

    /// find max 2 chains from term
    ch_score_t max_scores[2];
    ch_size_t bt_ids[2];
    max_scores[0] = max_scores[1] = CH_SCORE_MIN;
    printf("y, scores: ");
    for (ch_size_t i = 0; i < num_anchors; ++i) {
        ch_score_t score = scores[i + 1];
        printf("{%lu, %d}, ", anchors[i].y, score);
        if (score > max_scores[0]) {
            max_scores[0] = score;
            bt_ids[0] = i + 1;
        } else if (score > max_scores[1]) {
            max_scores[1] = score;
            bt_ids[1] = i + 1;
        }
    }
    printf("\n");
    printf("best scores: %d, %d\n", max_scores[0], max_scores[1]);
    printf("best ids: %d, %d\n", bt_ids[0], bt_ids[1]);
    //    printf("best ys: %lu, %lu\n", anchors[ids[0] - 1].y, anchors[ids[1] -
    //    1].y);

    /// back-tracking
    ch_result_t result;
    result.chains[0] = back_tracking(ids, bt_ids[0]);
    auto pos =
        anchors[result.chains[0].pos].x - anchors[result.chains[0].pos].y;
    result.chains[0].pos = pos;
    result.chains[0].score = max_scores[0];
    result.chains[1] = back_tracking(ids, bt_ids[1]);
    pos = anchors[result.chains[1].pos].x - anchors[result.chains[1].pos].y;
    result.chains[1].pos = pos;
    result.chains[1].score = max_scores[1];

    return result;
}
