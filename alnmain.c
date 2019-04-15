//
// Created by lisanhu on 9/26/18.
//

#include <zlib.h>
#include <time.h>
#include <stdio.h>

#include "accaln.h"
#include "kseq.h"
#include "mlog/logger.h"
#include "mutils.h"
#include "histo/histo.h"
#include "edlib/edlib.h"
#include "alnmain.h"


const double ERROR_RATE = 0.05;


static void _rev_comp_in_place(char *seq, uint32_t len);

void _rev_comp_in_place(char *seq, uint32_t len) {
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
        strncpy(buf + (max_len + 1) * i, reads[i].seq, reads[i].len);
        free(reads[i].seq);
        reads[i].seq = buf + (max_len + 1) * i;
    }
    return buf;
}

/// duplicate version in asindex.c
int32_t _dna_rand_ch() {
    static int32_t val = 0;
    static int pos = -1;
    if (pos < 0) {
        val = (int32_t) lrand48();
        pos = 0;
    } else if (pos < 31) {
        pos += 2;
    } else {
        val = (int32_t) lrand48();
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
            return i;
        }

        fread(&result[i].offset, sizeof(uint64_t), 1, mfp);
        fread(&result[i].seq_len, sizeof(size_t), 1, mfp);
    }
    return 0;
}


typedef struct seq_meta {
    uint64_t loc, off;
    mstring g_name;
    uint8_t strand;
} seq_meta;

#pragma acc routine seq
int
seq_lookup(const mta_entry *table, int len, uint64_t loc, uint32_t qlen,
           seq_meta *result) {
#pragma acc loop seq
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

    log.mvlog(&log, "fmi_read @ %s:%d", __FILE__, __LINE__);

    fmi_read(ctx->fmi, ctx->prefix);
    log.mvlog(&log, "fmi_read done.");

    char *path = cstr_concat(ctx->prefix, ".lch");
    log.mvlog(&log, "lc_read @ %s:%d", __FILE__, __LINE__);

    lc_read(path, ctx->lch);
    log.mvlog(&log, "lc_read done.");
    free(path);

    path = cstr_concat(ctx->genome, ".mta");
    mta_entry *meta = malloc(sizeof(mta_entry) * 65535);
    log.mvlog(&log, "ld_mta @ %s:%d", __FILE__, __LINE__);
    int len = load_mta(path, meta);
    log.mvlog(&log, "ld_mta done.");
    gen_sam_header(meta, len, stdout);
    ctx->mta = meta;
    ctx->mta_len = len;

    free(path);
    log.mvlog(&log, "ld_params @ %s:%d", __FILE__, __LINE__);
    params p = read_params("params", ctx);
    log.mvlog(&log, "ld_params done.");

    ctx->batch_size = p.batch_size; // default 1M reads
    ctx->histo_cap = p.thres;
    ctx->uninformative_thres = ctx->histo_cap;
    ctx->seed_len = p.seed_len;
    ctx->read2 = NULL;
//	ctx->sa_cache_sz = 1L << 29; // 32 / 8G x 8 Bytes
    ctx->sa_cache_sz = 1L << 33; // 8G entries
//	ctx->sa_cache_sz = 10000;

    u64 l;
    log.mvlog(&log, "load_file cat @ %s:%d", __FILE__, __LINE__);
    ctx->content = load_file(ctx->prefix, &l);
    log.mvlog(&log, "load_file done.");
    ctx->con_len = l;


    log.mvlog(&log, "Loading sa5 from disk. at %s:%d", __FILE__, __LINE__);
    char *fname = cstr_concat(ctx->prefix, ".sa5");
    FILE *stream = fopen(fname, "r");
    log.mvlog(&log, "Before ui40_read. at %s:%d", __FILE__, __LINE__);
    log.mvlog(&log, "%ld", l * sizeof(ui40_t));
    sa_buf = calloc(1, sizeof(sa_mem));
    sa_buf->mem = malloc(sizeof(ui40_t) * l);
    log.mvlog(&log, "Start ui40_read from disk. at %s:%d", __FILE__, __LINE__);
    sa_buf->len = ui40_fread(sa_buf->mem, l, stream);
    log.mvlog(&log, "Done ui40_read.");
    fclose(stream);
    free(fname);

    srand48(time(NULL));
}

#pragma acc routine seq
void remove_n(read_t *r) {
    const char *alpha = "ACGT";
#pragma acc loop seq
    for (u32 i = 0; i < r->len; ++i) {
        char ch = r->seq[i];
        if (ch == 'N' || ch == 'n') {
            r->seq[i] = alpha[_dna_rand_ch()];
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
    timer = log.mvlog(&log, "Done initializing, begin loading reference file %s",
              ctx.genome);

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
#pragma acc parallel loop
        for (u64 i = 0; i < len; ++i) {
            read_t r = reads[i];
            remove_n(&r);  /// todo: is this required?
            results[i].valid = false;


            histo *ot_iter_histo = histo_init(ctx.histo_cap);
            const int sl = ctx.seed_len; // todo: seed length should be further computed
            const int gl = 1;   // todo: gap length should be further computed

            entry best;
            double score = 0;
#pragma acc loop seq
            for (int iter = 0; iter < sl + gl; ++iter) {

                entry cand[2];
                histo *in_iter_histo = histo_init(ctx.histo_cap);
#pragma acc loop seq
                for (int j = iter; j < r.len - sl; j += sl + gl) {
                    u64 kk = 1, ll = ctx.fmi->length - 1, rr;

                    rr = lc_aln(r.seq + j, ctx.seed_len, &kk, &ll, ctx.fmi,
                                ctx.lch);

                    if (rr > 0 && rr < ctx.uninformative_thres) {
#pragma acc loop seq
                        for (u64 k = kk; k <= ll; ++k) {
                            u64 l = sa_access(ctx.prefix, ctx.sa_cache_sz, k) -
                                    j;
                            histo_add(in_iter_histo, l);
                        }
                    }
                }

                int num_seeds = r.len / (sl + gl);

                if (num_seeds > 0) {
                    u64 v = histo_find_2_max(in_iter_histo, cand);
                    score = (double) v / num_seeds;

//					double score = (double) v / num_seeds;
                    if (score >
                        0.6) { // todo: think of reasoning behind this threshold
                        // reason maybe the rest ratio are supposed to be around error rate
                        // todo: current result only support 1-1, need to think of other cases
                        best = cand[0];
                        break;
                    } else {
                        if (cand[0].val != 0) {
                            histo_add(ot_iter_histo, cand[0].key);
                        }
                    }
                }

                if (iter == sl + gl - 1) {
                    // last iteration
                    u64 v = histo_find_2_max(ot_iter_histo, cand);
                    best = cand[0];
                }

                histo_destroy(in_iter_histo);
            }


            u64 loc = best.key;
            int limit = (int) (ERROR_RATE * r.len * 2);
//            int limit = -1;
            seq_meta m;
            int meta_r = seq_lookup(ctx.mta, ctx.mta_len, loc, r.len, &m);
            if (m.strand == 1) {
                _rev_comp_in_place(r.seq, r.len);
            }


            char *cigar = cigar_align(r.seq, r.len, ctx.content + m.loc, r.len,
                                      &limit);
            result re = {.loc = loc, .off = m.off, .r_off = loc, .CIGAR = cigar,
                    .q_name = strdup(r.name.s), .g_name = strdup(
                            m.g_name.s), .qual = strdup(r.qual),
                    .query = strdup(
                            r.seq), .r_name = "*", .ed = limit, .mapq = 255,
                    .valid = (limit >= 0)};
//                    .valid = true};

            if (meta_r == 0) {
                re.valid = false;
            }
            results[i] = re;
            histo_destroy(ot_iter_histo);

        }

        log.mvlog(&log, "Done processing current batch, "
                        "currently processed %ld queries", total);


        FILE *out_stream = stdout;
        /// step 4: SAM generation
        for (int i = 0; i < len; ++i) {
            if (results[i].valid) {
                fprintf(out_stream,
                        "%s\t%d\t%s\t%ld\t%d\t%s\t%s\t%ld\t%d\t%s\t%s\tED:I:%d\n",
                        results[i].q_name, results[i].flag, results[i].g_name,
                        results[i].off + 1, results[i].mapq, results[i].CIGAR,
                        results[i].r_name, results[i].r_off, 0,
                        results[i].query,
                        results[i].qual, results[i].ed);
                free(results[i].CIGAR);
                valid += 1;
            }
        }
//		fclose(out_stream);
        free(buf);
        clock_gettime(CLOCK_MONOTONIC, &timer);
    }


    log.mvlog(&log, "Done aligning");
    log.mvlog(&log, "Sensitivity: %ld/%ld=%lf\n", valid, total,
              ((double) valid / total));


    kseq_destroy(seq);
    gzclose(fp);

    context_destroy(&ctx);

    return 0;
}


static inline int pair_end(int argc, const char *argv[]) {
    /// todo: need to implement this
    return -1;
}


int main(int argc, const char **argv) {

    if (argc != 3 && argc != 4) {
        usage(argv[0]);
        return -1;
    }

    if (argc == 3) {
        return single_end(argc, argv);
    }

    return pair_end(argc, argv);
}

params read_params(const char *path, context *ctx) {
    params result;
    FILE *fp = fopen(path, "r");
    result.thres = 300;
    result.batch_size = 1000000;
    result.seed_len = 20;
    if (fp) {
        fscanf(fp, "%lu %u %u", &result.batch_size, &result.seed_len,
               &result.thres);
    }

    mlog log = ctx->log;
    log.mvlog(&log, "Current settings:");
    log.mvlog(&log, "batch_size: %ld", result.batch_size);
    log.mvlog(&log, "seed_length: %d", result.seed_len);
    log.mvlog(&log, "non-informative seeds threshold: %d", result.thres);
    return result;
}
