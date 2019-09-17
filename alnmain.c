//
// Created by lisanhu on 9/26/18.
//

#include <stdio.h>
#include <time.h>
#include <zlib.h>

#include "accaln.h"
#include "alnmain.h"
#include "edlib/edlib.h"
#include "histo/histo.h"
#include "kseq.h"
#include "mlog/logger.h"
#include "mutils.h"

//#define CHUNK_SIZE 5000
#define CHUNK_SIZE 500
#define GE (-2)
#define EM (2)

const double ERROR_RATE = 0.05;
typedef uint64_t pos_t;

typedef struct _paired_end_pos_t {
    pos_t p1, p2;
} paired_end_pos_t;
typedef struct _anchor_t {
    pos_t x, y, w;
} _anchor_t;
_anchor_t const _anchor_zero = {0, 0, 0};
typedef struct _chain_resut {
    pos_t pos;
    int score;
} _chain_result;

typedef struct {
    u32 x, y;
} _2d_u32;

#pragma acc routine seq

static void _rev_comp_in_place(char *seq, uint32_t len);

void _rev_comp_in_place(char *seq, uint32_t len) {
    //#pragma acc loop seq
    for (uint32_t i = 0; i < len; ++i) {
        char c = seq[i];
        switch (c) {
            case 'A':
            case 'a':
                seq[i] = 'T';
                break;
            case 'C':
            case 'c':
                seq[i] = 'G';
                break;
            case 'G':
            case 'g':
                seq[i] = 'C';
                break;
            case 'T':
            case 't':
                seq[i] = 'A';
                break;
            default:
                /// should never come here
                seq[i] = 'N';
                break;
        }
    }

    for (uint32_t i = 0; i < len / 2; ++i) {
        char c = seq[i];
        seq[i] = seq[len - 1 - i];
        seq[len - 1 - i] = c;
    }
}

void gen_sam_header(mta_entry *mta, int l, FILE *stream) {
    long rg_id = time(NULL);
    char name[1024];

    for (int i = 0; i < l; ++i) {
        mta_entry m = mta[i];
        /// todo: assert m.seq_name.l < 1023
        strncpy(name, m.seq_name.s, m.seq_name.l);
        name[m.seq_name.l] = 0;
        fprintf(stream, "@SQ\tSN:%s\tLN:%ld\n", name, m.seq_len);
    }
    fprintf(stream, "@RG\tID:%s%ld\tSM:SM_data\n", "accaln", rg_id);
    fprintf(stream, "@PG\tID:%s\tPN:%s\n", "accaln", "accaln");
}

/**
 * A function to concat the storage of the reads together
 *
 * Note: After refactoring, reads could not be destroyed as normal
 * @param reads The reads array
 * @param l 	The number of reads in the array
 * @param ctx 	The program context to record max_len of the reads
 * @return 		The buffer used to store the reads
 */
char *refactor_reads_seq(read_t *reads, size_t l, context *ctx) {
    u32 max_len = 0;
    for (size_t i = 0; i < l; ++i) {
        max_len = max_len < reads[i].len ? reads[i].len : max_len;
    }
    ctx->max_read_len = max_len;

    char *buf = calloc((max_len + 1) * l, sizeof(char));

    for (size_t i = 0; i < l; ++i) {
        strncpy(buf + (max_len + 1) * i, reads[i].seq.s, reads[i].len);
        mstring_destroy(&reads[i].seq);
        reads[i].seq = mstring_from(buf + (max_len + 1) * i, false);
        //        reads[i].seq = buf + (max_len + 1) * i;
    }
    return buf;
}

/// duplicate version in asindex.c
int32_t _dna_rand_ch() {
    static int32_t val = 0;
    static int pos = -1;
    if (pos < 0) {
        val = (int32_t)lrand48();
        pos = 0;
    } else if (pos < 31) {
        pos += 2;
    } else {
        val = (int32_t)lrand48();
        pos = 0;
    }
    return (val >> pos) & 0x3;
}

/// todo: remove the assumption of less than 65535 genes
/// todo: sizeof(size_t) could be different on different systems
///	may result in indexing files could not be used in other systems
int load_mta(const char *path, mta_entry *result) {
    FILE *mfp = fopen(path, "r");
    for (int i = 0; i < 65535; ++i) {
        size_t l = mstring_read(&result[i].seq_name, mfp);
        if (l == 0) {
            //            mstring_destroy(&result[i].seq_name);
            fclose(mfp);
            return i;
        }

        fread(&result[i].offset, sizeof(uint64_t), 1, mfp);
        fread(&result[i].seq_len, sizeof(size_t), 1, mfp);
    }
    fclose(mfp);
    return 0;
}

typedef struct seq_meta {
    uint64_t loc, off;
    mstring g_name;
    uint8_t strand;
} seq_meta;

#pragma acc routine seq

int seq_lookup(const mta_entry *table, int len, uint64_t loc, uint32_t qlen,
               seq_meta *result) {
    //#pragma acc loop seq
    for (int i = 0; i < len; ++i) {
        uint64_t sl = table[i].seq_len;
        uint64_t start = table[i].offset;
        uint64_t end = start + sl * 2;
        if (loc >= start && loc + qlen <= start + sl) {
            /// the strand for genome
            result->strand = 0;
            result->g_name = table[i].seq_name;
            result->loc = loc;
            result->off = loc - start;
            return 1;
        } else if (loc >= start + sl && loc + qlen <= end) {
            /// the other strand for genome
            result->strand = 1;
            result->g_name = table[i].seq_name;
            result->off = end - loc - qlen;
            result->loc = result->off + start;
            return 1;
        }
    }
    return 0;
}

void init(context *ctx, int argc, const char **argv) {
    mlog log = ctx->log;

    ctx->genome = strdup(argv[1]);
    ctx->read1 = strdup(argv[2]);
    ctx->prefix = cstr_concat(ctx->genome, ".cat");

    ctx->fmi = malloc(sizeof(dna_fmi));
    ctx->lch = malloc(sizeof(lc_hash));

    fmi_read(ctx->fmi, ctx->prefix);
    log.mvlog(&log, "fmi_read done.");

    char *path = cstr_concat(ctx->prefix, ".lch");

    lc_read(path, ctx->lch);
    log.mvlog(&log, "lc_read done.");
    free(path);

    path = cstr_concat(ctx->genome, ".mta");
    mta_entry *meta = malloc(sizeof(mta_entry) * 65535);
    int len = load_mta(path, meta);
    log.mvlog(&log, "ld_mta done.");
    gen_sam_header(meta, len, stdout);
    ctx->mta = meta;
    ctx->mta_len = len;

    free(path);
    params p;
    if (argc != 6)
        p = read_params("params");
    else {
        p.batch_size = strtol(argv[3], NULL, 0);
        p.seed_len = strtol(argv[4], NULL, 0);
        p.thres = strtol(argv[5], NULL, 0);
    }
    log.mvlog(&log, "Current settings:");
    log.mvlog(&log, "batch_size: %ld", p.batch_size);
    log.mvlog(&log, "seed_length: %d", p.seed_len);
    log.mvlog(&log, "non-informative seeds threshold: %d", p.thres);
    log.mvlog(&log, "ld_params done.");

    ctx->batch_size = p.batch_size;  // default 1M reads
    ctx->histo_cap = p.thres;
    ctx->uninformative_thres = ctx->histo_cap;
    ctx->seed_len = p.seed_len;
    ctx->read2 = NULL;
    //	ctx->sa_cache_sz = 1L << 29; // 32 / 8G x 8 Bytes
    ctx->sa_cache_sz = 1L << 33;  // 8G entries
                                  //	ctx->sa_cache_sz = 10000;

    u64 l;
    ctx->content = load_file(ctx->prefix, &l);
    log.mvlog(&log, "load_file done.");
    ctx->con_len = l;

    char *fname = cstr_concat(ctx->prefix, ".sa5");
    FILE *stream = fopen(fname, "r");
    sa_buf = calloc(1, sizeof(sa_mem));
    sa_buf->mem = malloc(sizeof(ui40_t) * l);
    sa_buf->len = ui40_fread(sa_buf->mem, l, stream);
    log.mvlog(&log, "Done ui40_read.");
    fclose(stream);
    free(fname);

    srand48(time(NULL));
}

//#pragma acc routine seq
void remove_n(read_t *r) {
    const char *alpha = "ACGT";
    //#pragma acc loop seq
    for (u32 i = 0; i < r->len; ++i) {
        char ch = r->seq.s[i];
        if (ch == 'N' || ch == 'n') {
            r->seq.s[i] = alpha[_dna_rand_ch()];
        }
    }
}

static inline void usage(const char *path) {
    printf("Usage:\n");
    printf("\t%s ref.fa query.fq [query2.fq]\n", path);
}

static inline int single_end(int argc, const char *argv[]) {
    context ctx;
    ctx.log = new_mlogger(NULL);
    mlog log = ctx.log;
    struct timespec timer;

    log.mvlog(&log, "Start initialization");

    init(&ctx, argc, argv);
    timer = log.mvlog(
        &log, "Done initializing, begin loading reference file %s", ctx.genome);

    log.mvlog(&log, "Done loading reference in %lfs", time_elapse(timer));

    timer = log.mvlog(&log, "Begin loading queries from %s", ctx.read1);

    u64 len;
    gzFile fp = gzopen(ctx.read1, "r");
    kseq_t *seq = kseq_init(fp);
    read_t *reads = malloc(ctx.batch_size * sizeof(read_t));
    result *results = malloc(ctx.batch_size * sizeof(result));

    size_t total = 0, valid = 0;

    while ((len = reads_load(reads, ctx.batch_size, seq)) > 0) {
        char *buf = refactor_reads_seq(reads, len, &ctx);
        total += len;

        /// todo: loading part should use total length of bps instead of
        /// 		total number of reads
        log.mvlog(&log, "Done loading %ld queries in %lfs", len,
                  time_elapse(timer));

        log.mvlog(&log, "Begin processing queries");

        //        char **cigars = malloc(len * sizeof(char *));
        cigar *cig = malloc(len * sizeof(cigar));
        uint8_t **store = malloc(len * sizeof(uint8_t *));
        uint8_t *store_mem =
            malloc(len * ctx.max_read_len * 2 * sizeof(uint8_t));
        for (int i = 0; i < len; ++i) {
            store[i] = store_mem + i * ctx.max_read_len * 2;
        }

        for (int i = 0; i < len; ++i) {
            cig[i].cigar = store[i];
            cig[i].n_cigar_op = 0;
        }

        //#pragma acc parallel loop
        //#pragma omp parallel for
        entry best[CHUNK_SIZE];
        for (u64 i = 0; i < len; i += CHUNK_SIZE) {
            u64 max_limit = (i + CHUNK_SIZE > len) ? len - i : CHUNK_SIZE;

            for (u64 chunk_i = 0; chunk_i < max_limit; ++chunk_i) {
                read_t r = reads[i + chunk_i];
                //                remove_n(&r);  /// todo: is this required?
                results[i + chunk_i].valid = false;

                histo *ot_iter_histo = histo_init(ctx.histo_cap);
                const int sl = ctx.seed_len;  // todo: seed length should be
                                              // further computed
                const int gl =
                    1;  // todo: gap length should be further computed

                double score = 0;
                entry cand[2];
                int iter;

                for (iter = 0; iter < sl + gl; ++iter) {
                    histo *in_iter_histo = histo_init(ctx.histo_cap);

                    for (int j = iter; j < r.len - sl; j += sl + gl) {
                        u64 kk = 1, ll = ctx.fmi->length - 1, rr;

                        rr = lc_aln(r.seq.s + j, ctx.seed_len, &kk, &ll,
                                    ctx.fmi, ctx.lch);

                        if (rr > 0 && rr < ctx.uninformative_thres) {
                            for (u64 k = kk; k <= ll; ++k) {
                                u64 l =
                                    sa_access(ctx.prefix, ctx.sa_cache_sz, k) -
                                    j;
                                histo_add(in_iter_histo, l);
                            }
                        }
                    }

                    int num_seeds = r.len / (sl + gl);

                    if (num_seeds > 0) {
                        u64 v = histo_find_2_max(in_iter_histo, cand);
                        score = (double)v / num_seeds;

                        //					double score = (double) v /
                        // num_seeds;
                        if (score > 0.6) {  // todo: think of reasoning behind
                                            // this threshold
                            // reason maybe the rest ratio are supposed to be
                            // around error rate todo: current result only
                            // support 1-1, need to think of other cases
                            best[chunk_i] = cand[0];
                            histo_destroy(in_iter_histo);
                            break;
                        } else {
                            if (cand[0].val != 0) {
                                histo_add(ot_iter_histo, cand[0].key);
                            }
                        }
                    }

                    if (iter == sl + gl - 1) {
                        // last iteration
                        //                        u64 v =
                        //                        histo_find_2_max(ot_iter_histo,
                        //                        cand); best[chunk_i] =
                        //                        cand[0];
                    }

                    histo_destroy(in_iter_histo);
                }
                if (iter >= sl + gl - 1) {
                    histo_find_2_max(ot_iter_histo, cand);
                    best[chunk_i] = cand[0];
                }
                histo_destroy(ot_iter_histo);
            }

            ///// PART 2
            u64 loc[CHUNK_SIZE];
            seq_meta m[CHUNK_SIZE];
            int limit[CHUNK_SIZE];
            int meta_r[CHUNK_SIZE];

            mta_entry *mta = ctx.mta;
            int mta_len = ctx.mta_len;
            const char *content = ctx.content;
            int max_read_len = ctx.max_read_len;

            char *reads_mem = buf + (ctx.max_read_len + 1) * i;

#pragma acc parallel loop independent copyin(mta[:mta_len])            \
    copyin(content[:ctx.con_len])                                      \
        copy(reads_mem[:max_limit * (ctx.max_read_len + 1)])           \
            copy(limit[:], loc[:], meta_r[:], m[:], cig [i:max_limit]) \
                num_gangs(256) vector_length(256)
            for (u64 chunk_i = 0; chunk_i < max_limit; ++chunk_i) {
                read_t r = reads[i + chunk_i];
                loc[chunk_i] = best[chunk_i].key;
                limit[chunk_i] = (int)(ERROR_RATE * r.len * 2);
                //                int limit = -1;
                meta_r[chunk_i] =
                    seq_lookup(mta, mta_len, loc[chunk_i], r.len, m + chunk_i);
                if (m[chunk_i].strand == 1) {
                    //                        _rev_comp_in_place(
                    //                                reads_mem + (i + chunk_i)
                    //                                * (max_read_len + 1),
                    //                                r.len);
                    _rev_comp_in_place(r.seq.s, r.len);
                }

                cig[i + chunk_i] = cigar_align(
                    //                            reads_mem + (i + chunk_i) *
                    //                            (max_read_len + 1),
                    r.seq.s, r.len, content + m[chunk_i].loc, r.len,
                    &limit[chunk_i], store[chunk_i + i]);
                //                        store_mem + chunk_i * ctx.max_read_len
                //                        * 2);
                //                char cigar_buf[r.len * 2];
                //                parse_cigar(&cig[i + chunk_i], r.len,
                //                cigar_buf); printf("(i, chunk_i)=(%lu, %lu):
                //                %s\n", i, chunk_i, cigar_buf);
            }

            // PART 3
            /// todo: The query field may be different from original read
            /// because we use replace N in the reads and the mstring will
            /// update the original read data

            for (u64 chunk_i = 0; chunk_i < max_limit; ++chunk_i) {
                read_t r = reads[i + chunk_i];
                result re = {.loc = loc[chunk_i],
                             .off = m[chunk_i].off,
                             .r_off = loc[chunk_i],
                             .CIGAR = cig[i + chunk_i],
                             .q_name = r.name,
                             .g_name = m[chunk_i].g_name,
                             .qual = r.qual,
                             .query = r.seq,
                             .r_name = mstring_borrow("*", 1),
                             .ed = limit[chunk_i],
                             .mapq = 255,
                             .valid = (limit[chunk_i] >= 0),
                             .flag = 0};

                if (meta_r[chunk_i] == 0 || limit[chunk_i] == -1) {
                    re.valid = false;
                    re.flag += 0x4;
                    re.mapq = 0;
                } else {
                    if (m[chunk_i].strand == 1) {
                        re.flag += 16;
                    }
                }
                results[chunk_i + i] = re;
            }
        }

        log.mvlog(&log,
                  "Done processing current batch, "
                  "currently processed %ld queries",
                  total);

        FILE *out_stream = stdout;
        setvbuf(out_stream, NULL, _IOFBF, 4194304);  // 4MB buffer
        /// step 4: SAM generation
        for (int i = 0; i < len; ++i) {
            if (results[i].valid) {
                valid += 1;
            }

            if (results[i].query.l * 2 <= 0) {
                log.melog(&log, "Invalide read: %lu", results[i].query.l);
            }
            //            char *cigar_buf = malloc(results[i].query.l * 2);
            char cigar_buf[results[i].query.l * 2];
            parse_cigar(&results[i].CIGAR, results->query.l, cigar_buf);

            fprintf(out_stream,
                    "%.*s\t"  // query_name
                    "%d\t"    // flag
                    "%.*s\t"  // gene_name
                    "%ld\t"   //? results[i].off + 1
                    "%d\t"    // mapping quality
                    //                    "%.*s\t"        //CIGAR
                    "%s\t"        // CIGAR
                    "%.*s\t"      //??
                    "%ld\t"       // ?
                    "%d\t"        //?
                    "%.*s\t"      // query
                    "%.*s\t"      // quality string
                    "ED:I:%d\n",  // comment
                    (int)results[i].q_name.l, results[i].q_name.s,
                    results[i].flag, (int)results[i].g_name.l,
                    results[i].g_name.s, results[i].off + 1, results[i].mapq,
                    cigar_buf, (int)results[i].r_name.l, results[i].r_name.s,
                    0L, 0, (int)results[i].query.l, results[i].query.s,
                    (int)results[i].qual.l, results[i].qual.s, results[i].ed);
        }
        fflush(out_stream);
        free(cig);
        free(store);
        free(store_mem);
        reads_destroy(reads, len);
        free(buf);
        clock_gettime(CLOCK_MONOTONIC, &timer);
    }

    free(reads);
    free(results);

    log.mvlog(&log, "Done aligning");
    log.mvlog(&log, "Sensitivity: %ld/%ld=%lf\n", valid, total,
              ((double)valid / total));

    kseq_destroy(seq);
    gzclose(fp);

    context_destroy(&ctx);

    return 0;
}

int _chain_gap(_anchor_t a, _anchor_t b, int ge) {
    if (a.w == 0 || b.w == 0) {
        return 0;
    }
    return (b.x - a.x - a.w) * ge;
}

/**
 *  Report the id (seed_id, anchor_id) of best anchor which precedes anchor
 *  Because the gap scores are negative, we are looking for the smallest of its
 * complementary
 */
_2d_u32 max_gap_score_anchor(_anchor_t anchor, _anchor_t *anchors, u32 thres,
                             int seed_id, int *score) {
    _2d_u32 id = {0, 0};
    _anchor_t(*_anchors)[thres] = (_anchor_t(*)[thres])anchors;
    int min_rscore = 0;
    for (int sid = 0; sid < seed_id; ++sid) {
        for (int anchor_id = 0; anchor_id < thres; ++anchor_id) {
            int rscore = -_chain_gap(_anchors[sid][anchor_id], anchor, GE);
            if (rscore < min_rscore) {
                min_rscore = rscore;
                id.x = anchor_id;
                id.y = sid;
            }
        }
    }

    *score = -min_rscore;
    return id;
}

_chain_result chaining(_anchor_t *anchors, u32 thres, int num_seeds) {
    typedef struct {
        int score;
        _2d_u32 id;
    } cell;
    _chain_result result = {0, 0};
    if (num_seeds <= 0) return result;

    _anchor_t(*_anchors)[thres] = (_anchor_t(*)[thres])anchors;
    /// last 2 columns are used to record the pos of previous best anchor
    /// last row is used to store scores for terminate anchor
    cell scores[num_seeds + 1][thres + 1];
    for (u32 anchor_id = 0; anchor_id < thres; ++anchor_id) {
        scores[0][anchor_id].score = _anchors[0][anchor_id].w * EM;
        scores[0][anchor_id].id.x = 0;
        scores[0][anchor_id].id.y = 0;
    }

    for (int seed_id = 0; seed_id < num_seeds; ++seed_id) {
        for (u32 anchor_id = 0; anchor_id < thres; ++anchor_id) {
            scores[seed_id][anchor_id].score =
                _anchors[seed_id][anchor_id].w * EM;
            int gap_score;
            _2d_u32 id =
                max_gap_score_anchor(_anchors[seed_id][anchor_id], anchors,
                                     thres, seed_id, &gap_score);
            scores[seed_id][anchor_id].score += gap_score;
            scores[seed_id][anchor_id].id = id;
        }
    }
    int seed_id = num_seeds;
    int max_score = INT_MIN;
    _2d_u32 mid = {0, 0};
    for (u32 anchor_id = 0; anchor_id < thres; ++anchor_id) {
        scores[seed_id][anchor_id].score = 0;
        int gap_score;
        _anchor_t t = {0, 0, 0};
        _2d_u32 id =
            max_gap_score_anchor(t, anchors, thres, seed_id, &gap_score);
        scores[seed_id][anchor_id].score += gap_score;
        scores[seed_id][anchor_id].id = id;
        if (scores[seed_id][anchor_id].score > max_score) {
            max_score = scores[seed_id][anchor_id].score;
            mid = id;
        }
    }

    _2d_u32 best = mid;
    _2d_u32 last = best;

    while (scores[best.y][best.x].id.x != 0 &&
           scores[best.y][best.x].id.y != 0) {
        last = best;
        best = scores[best.y][best.x].id;
    }
    best = last;
    result.pos = _anchors[best.y][best.x].x;
    result.score = max_score;

    return result;
}

static paired_end_pos_t search_paired_reads(read_t *read1, read_t *read2,
                                            context *ctx) {
    paired_end_pos_t pos;
    histo *ot_iter_histo = histo_init(ctx->histo_cap);
    // todo: seed length should be further computed
    const int sl = ctx->seed_len;
    //    const int gl = 1;  // todo: gap length should be further computed

    double score = 0;
    _chain_result results[2];
    int iter;

    /// From minimap2 paper, the anchor is a triple (x, y, w)
    /// where x is the end pos in reference
    ///     y is the end pos in read
    ///     w is the seed length
    /// we will be using start pos for x and y
    for (iter = 0; iter < sl; ++iter) {
        //        histo *in_iter_histo = histo_init(ctx->histo_cap);
        int num_seeds = (read1->len - iter) / sl;
        u32 num_anchors = ctx->uninformative_thres * num_seeds;
        _anchor_t anchors[num_anchors];
        for (int i = 0; i < num_anchors; ++i) {
            anchors[i] = _anchor_zero;
        }

        for (int j = iter; j < read1->len - sl; j += sl) {
            u64 kk = 1, ll = ctx->fmi->length - 1, rr;
            rr = lc_aln(read1->seq.s + j, ctx->seed_len, &kk, &ll, ctx->fmi,
                        ctx->lch);
            if (rr > 0 && rr < ctx->uninformative_thres) {
                for (u64 k = kk; k <= ll; ++k) {
                    u64 l = sa_access(ctx->prefix, ctx->sa_cache_sz, k) - j;
                    int seed_id = (j - iter) / sl;
                    int ancho_id = seed_id * ctx->uninformative_thres + ll - kk;
                    anchors[ancho_id].x = l + j;
                    anchors[ancho_id].y = j;
                    /// seed length constant right now, use for future extension
                    anchors[ancho_id].w = sl;
                }
            }
        }

        results[0] = chaining(anchors, ctx->uninformative_thres, num_seeds);
    }
    for (iter = 0; iter < sl; ++iter) {
        //        histo *in_iter_histo = histo_init(ctx->histo_cap);
        int num_seeds = (read2->len - iter) / sl;
        u32 num_anchors = ctx->uninformative_thres * num_seeds;
        _anchor_t anchors[num_anchors];
        for (int i = 0; i < num_anchors; ++i) {
            anchors[i] = _anchor_zero;
        }

        for (int j = iter; j < read2->len - sl; j += sl) {
            u64 kk = 1, ll = ctx->fmi->length - 1, rr;
            rr = lc_aln(read2->seq.s + j, ctx->seed_len, &kk, &ll, ctx->fmi,
                        ctx->lch);
            if (rr > 0 && rr < ctx->uninformative_thres) {
                for (u64 k = kk; k <= ll; ++k) {
                    u64 l = sa_access(ctx->prefix, ctx->sa_cache_sz, k) - j;
                    int seed_id = (j - iter) / sl;
                    int ancho_id = seed_id * ctx->uninformative_thres + ll - kk;
                    anchors[ancho_id].x = l + j;
                    anchors[ancho_id].y = j;
                    /// seed length constant right now, use for future extension
                    anchors[ancho_id].w = sl;
                }
            }
        }

        results[1] = chaining(anchors, ctx->uninformative_thres, num_seeds);
    }

    /// todo: need to consider both reads in different helix
    /// need to consider a threshold when we choose to suppress one result and
    /// choose the other need to consider a secondary chain need to assign a
    /// score for insert size
    pos.p1 = results[0].pos;
    pos.p2 = results[1].pos;

    return pos;
}

static inline int pair_end(int argc, const char *argv[]) {
    /// todo: need to implement this
    context ctx;
    ctx.log = new_mlogger(NULL);
    mlog log = ctx.log;
    struct timespec timer;

    log.mvlog(&log, "Start initialization");

    init(&ctx, argc, argv);
    timer = log.mvlog(
        &log, "Done initializing, begin loading reference file %s", ctx.genome);

    log.mvlog(&log, "Done loading reference in %lfs", time_elapse(timer));

    timer = log.mvlog(&log, "Begin loading queries from %s", ctx.read1);

    u64 len;
    gzFile fp1 = gzopen(ctx.read1, "r");
    gzFile fp2 = gzopen(ctx.read2, "r");
    kseq_t *seq1 = kseq_init(fp1);
    kseq_t *seq2 = kseq_init(fp1);
    read_t *reads1 = malloc(ctx.batch_size * sizeof(read_t));
    read_t *reads2 = malloc(ctx.batch_size * sizeof(read_t));
    result *result1 = malloc(ctx.batch_size * sizeof(result));
    result *result2 = malloc(ctx.batch_size * sizeof(result));

    size_t total = 0, valid = 0, batch_size = ctx.batch_size;

    while ((len = reads_load(reads1, batch_size, seq1)) ==
           reads_load(reads2, batch_size, seq2)) {
        char *buf1 = refactor_reads_seq(reads1, len, &ctx);
        char *buf2 = refactor_reads_seq(reads2, len, &ctx);
        total += len;
        log.mvlog(&log, "Done loading %ld paired reads in %lfs", len,
                  time_elapse(timer));

        /// allocate memory for cigar strings
        /// use stack memory because easier to port to GPU
        /// the use of stack memory restrict the batch_size,
        /// but should be fine
        cigar cig1[len];
        cigar cig2[len];
        uint8_t store1[len][len * ctx.max_read_len * 2];
        uint8_t store2[len][len * ctx.max_read_len * 2];
        for (int i = 0; i < len; ++i) {
            cig1[i].cigar = store1[i];
            cig2[i].cigar = store2[i];
            cig1[i].n_cigar_op = 0;
            cig2[i].n_cigar_op = 0;
        }

        //#pragma acc parallel loop
        //#pragma omp parallel for
        entry best[batch_size];

        mta_entry *mta = ctx.mta;
        int mta_len = ctx.mta_len;
        for (u64 i = 0; i < len; ++i) {
            read_t *r1 = &reads1[i];
            read_t *r2 = &reads2[i];

            paired_end_pos_t p = search_paired_reads(r1, r2, &ctx);
            int m[2];
            m[0] = seq_lookup(mta,mta_len, p.p1, r1->len);
            m[1] = seq_lookup(mta,mta_len, p.p2, r2->len);
        }

        free(buf1);
        free(buf2);
    }
    //
    //        for (u64 i = 0; i < len; i += CHUNK_SIZE) {
    //            u64 max_limit = (i + CHUNK_SIZE > len) ? len - i : CHUNK_SIZE;
    //
    //
    //            ///// PART 2
    //            u64 loc[CHUNK_SIZE];
    //            seq_meta m[CHUNK_SIZE];
    //            int limit[CHUNK_SIZE];
    //            int meta_r[CHUNK_SIZE];
    //
    //            mta_entry *mta = ctx.mta;
    //            int mta_len = ctx.mta_len;
    //            const char *content = ctx.content;
    //            int max_read_len = ctx.max_read_len;
    //
    //            char *reads_mem = buf + (ctx.max_read_len + 1) * i;
    //
    //#pragma acc parallel loop independent copyin(mta[:mta_len])            \
//    copyin(content[:ctx.con_len])                                      \
//        copy(reads_mem[:max_limit * (ctx.max_read_len + 1)])           \
//            copy(limit[:], loc[:], meta_r[:], m[:], cig [i:max_limit]) \
//                num_gangs(256) vector_length(256)
    //            for (u64 chunk_i = 0; chunk_i < max_limit; ++chunk_i) {
    //                read_t r = reads[i + chunk_i];
    //                loc[chunk_i] = best[chunk_i].key;
    //                limit[chunk_i] = (int)(ERROR_RATE * r.len * 2);
    //                //                int limit = -1;
    //                meta_r[chunk_i] =
    //                    seq_lookup(mta, mta_len, loc[chunk_i], r.len, m
    //                    chunk_i);
    //                if (m[chunk_i].strand == 1) {
    //                    //                        _rev_comp_in_place(
    //                    //                                reads_mem + (i
    //                    chunk_i)
    //                    //                                * (max_read_len
    //                    1),
    //                    //                                r.len);
    //                    _rev_comp_in_place(r.seq.s, r.len);
    //                }
    //
    //                cig[i + chunk_i] = cigar_align(
    //                    //                            reads_mem + (i
    //                    chunk_i) *
    //                    //                            (max_read_len + 1),
    //                    r.seq.s, r.len, content + m[chunk_i].loc, r.len,
    //                    &limit[chunk_i], store[chunk_i + i]);
    //                //                        store_mem + chunk_i *
    //                ctx.max_read_len
    //                //                        * 2);
    //                //                char cigar_buf[r.len * 2];
    //                //                parse_cigar(&cig[i + chunk_i], r.len,
    //                //                cigar_buf); printf("(i, chunk_i)=(%lu,
    //                %lu):
    //                //                %s\n", i, chunk_i, cigar_buf);
    //            }
    //
    //            // PART 3
    //            /// todo: The query field may be different from original read
    //            /// because we use replace N in the reads and the mstring will
    //            /// update the original read data
    //
    //            for (u64 chunk_i = 0; chunk_i < max_limit; ++chunk_i) {
    //                read_t r = reads[i + chunk_i];
    //                result re = {.loc = loc[chunk_i],
    //                             .off = m[chunk_i].off,
    //                             .r_off = loc[chunk_i],
    //                             .CIGAR = cig[i + chunk_i],
    //                             .q_name = r.name,
    //                             .g_name = m[chunk_i].g_name,
    //                             .qual = r.qual,
    //                             .query = r.seq,
    //                             .r_name = mstring_borrow("*", 1),
    //                             .ed = limit[chunk_i],
    //                             .mapq = 255,
    //                             .valid = (limit[chunk_i] >= 0),
    //                             .flag = 0};
    //
    //                if (meta_r[chunk_i] == 0 || limit[chunk_i] == -1) {
    //                    re.valid = false;
    //                    re.flag += 0x4;
    //                    re.mapq = 0;
    //                } else {
    //                    if (m[chunk_i].strand == 1) {
    //                        re.flag += 16;
    //                    }
    //                }
    //                results[chunk_i + i] = re;
    //            }
    //        }
    //
    //        log.mvlog(&log,
    //                  "Done processing current batch, "
    //                  "currently processed %ld queries",
    //                  total);
    //
    //        FILE *out_stream = stdout;
    //        setvbuf(out_stream, NULL, _IOFBF, 4194304);  // 4MB buffer
    //        /// step 4: SAM generation
    //        for (int i = 0; i < len; ++i) {
    //            if (results[i].valid) {
    //                valid += 1;
    //            }
    //
    //            if (results[i].query.l * 2 <= 0) {
    //                log.melog(&log, "Invalide read: %lu", results[i].query.l);
    //            }
    //            //            char *cigar_buf = malloc(results[i].query.l *
    //            2); char cigar_buf[results[i].query.l * 2];
    //            parse_cigar(&results[i].CIGAR, results->query.l, cigar_buf);
    //
    //            fprintf(out_stream,
    //                    "%.*s\t"  // query_name
    //                    "%d\t"    // flag
    //                    "%.*s\t"  // gene_name
    //                    "%ld\t"   //? results[i].off + 1
    //                    "%d\t"    // mapping quality
    //                    //                    "%.*s\t"        //CIGAR
    //                    "%s\t"        // CIGAR
    //                    "%.*s\t"      //??
    //                    "%ld\t"       // ?
    //                    "%d\t"        //?
    //                    "%.*s\t"      // query
    //                    "%.*s\t"      // quality string
    //                    "ED:I:%d\n",  // comment
    //                    (int)results[i].q_name.l, results[i].q_name.s,
    //                    results[i].flag, (int)results[i].g_name.l,
    //                    results[i].g_name.s, results[i].off + 1,
    //                    results[i].mapq, cigar_buf, (int)results[i].r_name.l,
    //                    results[i].r_name.s, 0L, 0, (int)results[i].query.l,
    //                    results[i].query.s, (int)results[i].qual.l,
    //                    results[i].qual.s, results[i].ed);
    //        }
    //        fflush(out_stream);
    //        free(cig);
    //        free(store);
    //        free(store_mem);
    //        reads_destroy(reads, len);
    //        free(buf);
    //        clock_gettime(CLOCK_MONOTONIC, &timer);
    //    }

    free(reads1);
    free(reads2);
    free(result1);
    free(result2);

    log.mvlog(&log, "Done aligning");
    log.mvlog(&log, "Sensitivity: %ld/%ld=%lf\n", valid, total,
              ((double)valid / total));

    kseq_destroy(seq1);
    kseq_destroy(seq2);
    gzclose(fp1);
    gzclose(fp2);

    context_destroy(&ctx);

    return 0;
}

int main(int argc, const char **argv) {
    if (argc != 6 && argc != 3 && argc != 4) {
        usage(argv[0]);
        return -1;
    }

    if (argc == 3 || argc == 6) {
        return single_end(argc, argv);
    }

    return pair_end(argc, argv);
}

params read_params(const char *path) {
    params result;
    FILE *fp = fopen(path, "r");
    result.thres = 300;
    //    result.batch_size = 1000000;
    result.batch_size = 1000;
    result.seed_len = 20;
    if (fp) {
        fscanf(fp, "%lu %u %u", &result.batch_size, &result.seed_len,
               &result.thres);
        fclose(fp);
    }

    return result;
}
