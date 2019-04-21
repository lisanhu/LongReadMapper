//
// Created by lisanhu on 9/26/18.
//

#ifndef ACCSEQV8_ACCALN_H
#define ACCSEQV8_ACCALN_H


#include <stdint.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdbool.h>

#include "mutils.h"
#include "kseq.h"
#include "fmidx/fmidx.h"
#include "lchash/lchash.h"
#include "mlog/logger.h"

KSEQ_INIT(gzFile, gzread)


#ifdef __cplusplus
extern "C" {
#endif

typedef uint64_t u64;
typedef uint32_t u32;
//typedef int64_t i64;
//typedef int32_t i32;

#ifndef MFREE
#define MFREE(x) {if ((x)) {free((x)); (x) = NULL;}}
#endif


/**
* Definition of class read_t
*/
typedef struct {
	mstring seq, qual, name;
	size_t rid;
	u32 len;
} read_t;

/**
 * Destructor for read_t
 * @param rd	the pointer to the read_t object to be destroyed
 */
void read_destroy(read_t *rd);
/**
 * Batch destructor for read_t
 * @param rd 	the reads array to be cleared
 * @param num 	the length of the reads array
 */
void reads_destroy(read_t *rd, size_t num);
/**
 * Batch constructor for read_t
 * @param buf	the buffer to store the read_t objects
 * @param num	the maximum number of reads to parse
 * @param seq	the kseq_t buffer used to continuing reading the same file
 * @return		the number of reads actually consumed
 */
size_t reads_load(read_t *buf, size_t num, kseq_t *seq);

typedef struct mta_entry {
	mstring seq_name;
	uint64_t offset;
	size_t seq_len;
} mta_entry;

typedef struct _context {
	char *read1, *read2, *genome, *prefix;
	u64 batch_size, uninformative_thres, sa_cache_sz;
	u32 max_read_len, seed_len, histo_cap;
	dna_fmi *fmi;
	lc_hash *lch;
	const char *content;
	u64 con_len;
	mta_entry *mta;
	int mta_len;
	mlog log;
//	size_t batch_size, gene_total_len;
//	int32_t max_len, seed_len;
//	char *read1, *read2, *genome, *content, *out_file; /// content: concat file content
//	bwt_t *bwt;
//	bntseq_t *bns;
//	gnl_vec_t gnl_vec;
//	fmidx_dna *dna;
} context;

void context_destroy(context *ctx);

typedef struct result{
	/// location in the file, offset in BAM format, offset of mate/pair read
	uint64_t loc, off, r_off;
	/// CIGAR string, query name, gene name, query quality, query, reference name of mate/pair read
	mstring CIGAR, q_name, g_name, qual, query, r_name;
	int ed, mapq, flag;
	bool valid;
} result;


#ifdef __cplusplus
};
#endif

#endif //ACCSEQV8_ACCALN_H
