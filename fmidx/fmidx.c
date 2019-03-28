//
// Created by lisanhu on 9/17/18.
//

#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <stddef.h>
#include <unistd.h>  /// dup

#include "fmidx.h"
#include "../psascan/sa_use.h"
#include "../mutils.h"


sa_mem *sa_buf = NULL;


inline uint64_t sa_access(const char *prefix, uint64_t cache_sz, uint64_t loc) {
//	if (sa_buf && loc < sa_buf->start + sa_buf->len && loc >= sa_buf->start) {
//	} else {
//		/// when we need to load items from disk
//        printf("Loading from disk. on line %d\n", __LINE__);
//		if (!sa_buf) {
//			sa_buf = calloc(1, sizeof(sa_mem));
//			sa_buf->mem = malloc(cache_sz * sizeof(ui40_t));
//		}
//		/// initial sa_buf mem
//        printf("Before memset on line %d\n", __LINE__);
////		memset(sa_buf->mem, 0, cache_sz * sizeof(ui40_t));
//        printf("Done memset on line %d\n", __LINE__);
//		/// load the block for loc
//		uint64_t skipped = loc / cache_sz;
//		char *fname = cstr_concat(prefix, ".sa5");
//		FILE *stream = fopen(fname, "r");
//
//		fseek(stream, (long) (skipped * cache_sz * 5), SEEK_SET);
//        printf("Before ui40_fread at line %d\n", __LINE__);
//		sa_buf->len = ui40_fread(sa_buf->mem, cache_sz, stream);
//		sa_buf->start = skipped * cache_sz;
//        printf("Done ui40_fread at line %d\n", __LINE__);
//
//		free(fname);
//		fclose(stream);
//        printf("Done loading from disk. on line %d\n", __LINE__);
//	}
//	return ui40_convert(sa_buf->mem[loc - sa_buf->start]);
    if (sa_buf != NULL) {
        if (loc > sa_buf->len) {
            printf("Accessing sa element out of scope");
			return 0;
        }
    } else {
		sa_buf = calloc(1, sizeof(sa_mem));
		sa_buf->mem = malloc(sizeof(ui40_t) * cache_sz);
		char *fname = cstr_concat(prefix, ".sa5");
		FILE *stream = fopen(fname, "r");
		sa_buf->len = ui40_fread(sa_buf->mem, cache_sz, stream);
		free(fname);
    }
//    else {
//        printf("Loading sa5 from disk. at %s:%d\n", __FILE__, __LINE__);
//        char *fname = cstr_concat(prefix, ".sa5");
//        FILE *stream = fopen(fname, "r");
//        printf("Start ui40_read from disk. at %s:%d\n", __FILE__, __LINE__);
//        sa_buf->len = ui40_fread(sa_buf->mem, cache_sz, stream);
//        printf("Done ui40_read.\n");
//    }
    return ui40_convert(sa_buf->mem[loc]);
}

void sa_done_access() {
	free(sa_buf->mem);
	sa_buf = NULL;
}

char text_buf_access(uint64_t *start, uint64_t size, const char *prefix,
                     char *buf, uint64_t i, bool init) {
	uint64_t blkid = i / size;
	uint64_t offset = i % size;
	size_t items_read;

	if (init || *start / size != blkid) {
		FILE *fp = fopen(prefix, "r");

		if (!fp) {
			fprintf(stderr, "Cannot open file %s\n", prefix);
		}

		fseek(fp, (long) (blkid * size * sizeof(char)), SEEK_SET);
		items_read = fread(buf, sizeof(char), size, fp);
		*start = blkid * size;

		/// need to further consider index ot of bounds case
		if (items_read <= offset) {
			/// index out of bounds
		}
		fclose(fp);
	}
//
	return buf[offset];
}


/**
 * Build bwt string from sa5 file. Memory consumption: n + ram_use
 * n is the size of the prefix file
 * ram_use is the available buffer budget
 * @param prefix
 * @param ram_use
 * @return mstring format of the bwt string
 */
static mstring _bwt_from_sa5(const char *prefix, long ram_use) {
	mstring bwt;
	uint64_t size = (uint64_t) ram_use;

	char *buf = malloc(size * sizeof(char));
	bwt.l = file_length(prefix);
	bwt.s = malloc((bwt.l + 1) * sizeof(char));
	bwt.s[bwt.l] = '\0';

	uint64_t start = 0;
	text_buf_access(&start, size, prefix, buf, 0, true);

	for (uint64_t i = 0; i < bwt.l; ++i) {
		uint64_t si = sa_access(prefix, bwt.l, i);
		bwt.s[i] = (char) (si == 0 ? '$' :
		                   text_buf_access(&start, size, prefix, buf, si - 1,
		                                   false));
	}

	free(buf);

	return bwt;
}


static void _build_c_table(const char *prefix, uint64_t *tab, long ram_use) {
	size_t ram = (size_t) ram_use;
	uint64_t start = 0;
	char *buf = malloc(ram * sizeof(char));
	text_buf_access(&start, ram, prefix, buf, 0, true);
	size_t length = file_length(prefix);

	memset(tab, 0, 256 * sizeof(uint64_t));

	for (uint64_t i = 0; i < length - 1; ++i) {
		char ch = text_buf_access(&start, ram, prefix, buf, i, false);
//		if (0x41 <= ch && ch <= 0x5A) {
//			ch += 0x20;
//		}
		tab[ch]++;
	}

	uint64_t sum = 0;
	for (int i = 0; i < 256; ++i) {
		uint64_t tmp = sum + tab[i];
		tab[i] = sum;
		sum = tmp;
	}
	free(buf);
}


static void _build_o_dna_table(const char *bwt, uint64_t length, int ratio, uint64_t *tab) {
	int mapper[256];
	mapper['a'] = mapper['A'] = 0;
	mapper['c'] = mapper['C'] = 1;
	mapper['g'] = mapper['G'] = 2;
	mapper['t'] = mapper['T'] = 3;

	uint64_t tmp[4];
	tmp[0] = tmp[1] = tmp[2] = tmp[3] = 0;
	for (uint64_t i = 0; i < length; ++i) {
		uint64_t idx = i / ratio;
		uint64_t mod = i % ratio;

		if (mod == 0) {
			tab[4 * idx + 0] = tmp[0];
			tab[4 * idx + 1] = tmp[1];
			tab[4 * idx + 2] = tmp[2];
			tab[4 * idx + 3] = tmp[3];
		}
		char c = bwt[i];
		if (c != '\0' && c != '$') tmp[mapper[c]]++;
	}
}


void fmi_build(const char *prefix, dna_fmi *idx, int o_ratio, long ram_use) {
	/// build .sa5
	sa_build(prefix, ram_use);
	printf("Done building suffix array\n");

	/// build c
	printf("Start building C table\n");
	idx->c = malloc(256 * sizeof(uint64_t));
	if (idx->c == NULL) {
		printf("Fail to malloc 256 uint64_t\n");
	}
	_build_c_table(prefix, idx->c, ram_use);

	/// build bwt
	printf("Start building bwt from suffix array\n");
	mstring bwt = _bwt_from_sa5(prefix, ram_use);
	idx->bwt = bwt.s;
	idx->length = bwt.l;

	/// build o assuming only ACGT are found in the text
	printf("Start building O table\n");
	idx->o_ratio = o_ratio;
	idx->o_len = 4 * (idx->length / o_ratio + 1);
	idx->o = malloc(idx->o_len * sizeof(uint64_t));
	_build_o_dna_table(idx->bwt, idx->length, o_ratio, idx->o);
}


void fmi_destroy(dna_fmi *idx) {
	if (idx->c) {
		free(idx->c);
		idx->c = NULL;
	}

	if (idx->o) {
		free(idx->o);
		idx->o = NULL;
	}

	if (idx->bwt) {
		free(idx->bwt);
		idx->bwt = NULL;
	}

	idx->length = 0;
}


void fmi_write(dna_fmi idx, const char *prefix) {
	char *fname = cstr_concat(prefix, ".mfi");
	FILE *fp = fopen(fname, "w");

	/// write c
	fwrite(idx.c, sizeof(uint64_t), 256, fp);

	/// write o_ratio, o_len, o
	fwrite(&idx.o_ratio, sizeof(int), 1, fp);
	fwrite(&idx.o_len, sizeof(uint64_t), 1, fp);
	fwrite(idx.o, sizeof(uint64_t), idx.o_len, fp);

	/// write length, bwt
	fwrite(&idx.length, sizeof(uint64_t), 1, fp);
	fwrite(idx.bwt, sizeof(char), idx.length, fp);

	fclose(fp);
	free(fname);
}

void fmi_read(dna_fmi *idx, const char *prefix) {
	char *fname = cstr_concat(prefix, ".mfi");
	FILE *fp = fopen(fname, "r");

	/// read c
	idx->c = malloc(256 * sizeof(uint64_t));
	fread(idx->c, sizeof(uint64_t), 256, fp);

	/// read o_ratio, o_len, o
	fread(&idx->o_ratio, sizeof(int), 1, fp);
	fread(&idx->o_len, sizeof(uint64_t), 1, fp);
	idx->o = malloc(idx->o_len * sizeof(uint64_t));
	fread(idx->o, sizeof(uint64_t), idx->o_len, fp);

	/// read length, bwt
	fread(&idx->length, sizeof(uint64_t), 1, fp);
	idx->bwt = malloc(idx->length + 1);
	fread(idx->bwt, sizeof(char), idx->length, fp);
	idx->bwt[idx->length] = 0;

	fclose(fp);
	free(fname);
}

uint64_t _occ_access(const dna_fmi *idx, char c, uint64_t loc) {
	int mapper[256];
	mapper['a'] = mapper['A'] = 0;
	mapper['c'] = mapper['C'] = 1;
	mapper['g'] = mapper['G'] = 2;
	mapper['t'] = mapper['T'] = 3;

	int ratio = idx->o_ratio;
	// todo: move mapper to global
	uint64_t id = loc / ratio;
//	uint64_t mod = loc % ratio;
	uint64_t count = 0;
	for (uint64_t i = id * ratio; i <= loc; ++i) {
		if (idx->bwt[i] == c) count++;
	}
	return idx->o[4 * id + mapper[c]] + count;
}

uint64_t fmi_aln(const dna_fmi *idx, const char *qry, int len, uint64_t *k,
                 uint64_t *l) {
	int mapper[256];
	mapper['a'] = mapper['A'] = 0;
	mapper['c'] = mapper['C'] = 1;
	mapper['g'] = mapper['G'] = 2;
	mapper['t'] = mapper['T'] = 3;
	uint64_t kk = *k, ll = *l;
	for (int i = len - 1; i >= 0; --i) {
		char c = qry[i];
		// todo: move 128 to config
		kk = idx->c[c] + _occ_access(idx, c, kk - 1) + 1;
		ll = idx->c[c] + _occ_access(idx, c, ll);
		if (kk > ll) break;
	}
	*k = kk;
	*l = ll;
	return kk > ll ? 0 : ll - kk + 1;
}

