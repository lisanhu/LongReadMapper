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

//#define CHUNK_SIZE 5000
#define CHUNK_SIZE 500

const double ERROR_RATE = 0.05;

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

int
seq_lookup(const mta_entry *table, int len, uint64_t loc, uint32_t qlen,
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
    params p;
    if (argc != 6) p = read_params("params");
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
    timer = log.mvlog(&log,
                      "Done initializing, begin loading reference file %s",
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

//        char **cigars = malloc(len * sizeof(char *));
        cigar *cig = malloc(len * sizeof(cigar));
        uint8_t **store = malloc(len * sizeof(uint8_t *));
        uint8_t *store_mem = malloc(
                len * ctx.max_read_len * 2 * sizeof(uint8_t));
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
                const int sl = ctx.seed_len; // todo: seed length should be further computed
                const int gl = 1;   // todo: gap length should be further computed

                double score = 0;
                entry cand[2];
                int iter;


                for (iter = 0; iter < sl + gl; ++iter) {

                    histo *in_iter_histo = histo_init(ctx.histo_cap);

                    for (int j = iter; j < r.len - sl; j += sl + gl) {
                        u64 kk = 1, ll = ctx.fmi->length - 1, rr;

                        rr = lc_aln(r.seq.s + j, ctx.seed_len, &kk, &ll,
                                    ctx.fmi,
                                    ctx.lch);

                        if (rr > 0 && rr < ctx.uninformative_thres) {

                            for (u64 k = kk; k <= ll; ++k) {
                                u64 l = sa_access(ctx.prefix, ctx.sa_cache_sz,
                                                  k) -
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
                        //                        u64 v = histo_find_2_max(ot_iter_histo, cand);
                        //                        best[chunk_i] = cand[0];
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

#pragma acc parallel loop independent copyin(mta[:mta_len]) \
            copyin(content[:ctx.con_len]) \
            copy(reads_mem[:max_limit * (ctx.max_read_len + 1)]) \
            copy(limit[:], loc[:], meta_r[:], m[:], cig[i:max_limit]) \
                num_gangs(256) vector_length(256)
            for (u64 chunk_i = 0; chunk_i < max_limit; ++chunk_i) {
                read_t r = reads[i + chunk_i];
                loc[chunk_i] = best[chunk_i].key;
                limit[chunk_i] = (int) (ERROR_RATE * r.len * 2);
                //                int limit = -1;
                meta_r[chunk_i] = seq_lookup(mta, mta_len, loc[chunk_i],
                                             r.len,
                                             m + chunk_i);
                if (m[chunk_i].strand == 1) {
//                        _rev_comp_in_place(
//                                reads_mem + (i + chunk_i) * (max_read_len + 1),
//                                r.len);
                    _rev_comp_in_place(r.seq.s, r.len);
                }

                cig[i + chunk_i] = cigar_align(
//                            reads_mem + (i + chunk_i) * (max_read_len + 1),
                        r.seq.s,
                        r.len,
                        content + m[chunk_i].loc,
                        r.len, &limit[chunk_i],
                        store[chunk_i + i]);
//                        store_mem + chunk_i * ctx.max_read_len * 2);
//                char cigar_buf[r.len * 2];
//                parse_cigar(&cig[i + chunk_i], r.len, cigar_buf);
//                printf("(i, chunk_i)=(%lu, %lu): %s\n", i, chunk_i, cigar_buf);
            }

            // PART 3
            /// todo: The query field may be different from original read
            /// because we use replace N in the reads and the mstring will
            /// update the original read data

            for (u64 chunk_i = 0; chunk_i < max_limit; ++chunk_i) {
                read_t r = reads[i + chunk_i];
                result re = {.loc = loc[chunk_i], .off = m[chunk_i].off, .r_off = loc[chunk_i],
                        .CIGAR = cig[i + chunk_i], .q_name = r.name,
                        .g_name = m[chunk_i].g_name, .qual = r.qual, .query = r.seq,
                        .r_name = mstring_borrow("*", 1), .ed = limit[chunk_i],
                        .mapq = 255, .valid = (limit[chunk_i] >= 0), .flag = 0};

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

        log.mvlog(&log, "Done processing current batch, "
                        "currently processed %ld queries", total);


        FILE *out_stream = stdout;
        setvbuf(out_stream, NULL, _IOFBF, 4194304); // 4MB buffer
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
                    "%.*s\t"        //query_name
                    "%d\t"          //flag
                    "%.*s\t"        //gene_name
                    "%ld\t"         //? results[i].off + 1
                    "%d\t"          //mapping quality
                    //                    "%.*s\t"        //CIGAR
                    "%s\t"        //CIGAR
                    "%.*s\t"        //??
                    "%ld\t"         // ?
                    "%d\t"          //?
                    "%.*s\t"        //query
                    "%.*s\t"        //quality string
                    "ED:I:%d\n",    //comment
                    (int) results[i].q_name.l, results[i].q_name.s,
                    results[i].flag,
                    (int) results[i].g_name.l, results[i].g_name.s,
                    results[i].off + 1,
                    results[i].mapq,
                    cigar_buf,
                    (int) results[i].r_name.l, results[i].r_name.s,
                    0L,
                    0,
                    (int) results[i].query.l, results[i].query.s,
                    (int) results[i].qual.l, results[i].qual.s,
                    results[i].ed);
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
