//
// Created by lisanhu on 9/24/18.
//

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "lchash.h"
#include "../accaln.h"

inline uint64_t lc_access(lc_hash h, uint64_t loc, uint64_t *k, uint64_t *l) {
	*k = h.lc[2 * loc];
	*l = h.lc[2 * loc + 1];
	return *l - *k;
}


static const char *_seq_from_num(uint64_t num, int hlen) {
	if (hlen <= 0) return "";
	char buf[hlen + 1];
	buf[hlen] = '\0';
	const char *alphab = "ACGT";

	uint64_t mask = 3;

	for (int i = 0; i < hlen; ++i) {
		uint64_t tmp = num >> (2 * i) & mask;
		buf[hlen - 1 - i] = alphab[tmp];
	}

	return strdup(buf);
}


inline static uint64_t _num_from_seq(const char *seq, int hlen) {
	uint64_t sum = 0;
	int mapper[256];
	mapper['a'] = mapper['A'] = 0;
	mapper['c'] = mapper['C'] = 1;
	mapper['g'] = mapper['G'] = 2;
	mapper['t'] = mapper['T'] = 3;
	for (int i = 0; i < hlen; ++i) {
		sum += mapper[seq[i]];
		sum <<= 2;
		// todo: warning possible overflow
	}
	return sum >> 2;
}


lc_hash lc_build(const dna_fmi *fmi, int hlen) {
	uint64_t upper = lc_length(hlen);
	uint64_t *lc = malloc(sizeof(uint64_t) * upper * 2);
	for (uint64_t i = 0; i < upper; ++i) {
		uint64_t k = 1, l = fmi->length - 1;
		const char *seq = _seq_from_num(i, hlen);
		uint64_t r = fmi_aln(fmi, seq, hlen, &k, &l);
//		if (i == 10243342) {
//			printf("%s\n%lu\n%lu %lu\n", seq, r, k, l);
//		}
		free((void *) seq);

		if (r == 0) {
			k = l = 0;
		}
		lc[2 * i] = k;
		lc[2 * i + 1] = l;
	}

	lc_hash hash = {.lc = lc, .hlen = hlen, .len = 2 * upper};
	return hash;
}

uint64_t lc_length(int hlen) {
	return 1U << (2 * hlen);
}

void lc_destroy(lc_hash hash) {
	if (hash.lc == NULL) {
		return;
	}

	free(hash.lc);
//	lc_hash h = {.lc = NULL, .len = 0, .hlen = 0};
//	hash = h;
}

#pragma acc routine seq
inline uint64_t
lc_aln(const char *qry, int qlen, uint64_t *k, uint64_t *l, const dna_fmi *fmi,
       const lc_hash *hash) {

	int left = qlen - hash->hlen;
	if (qlen >= hash->hlen) {
		uint64_t num = _num_from_seq(qry + left, hash->hlen);
		lc_access(*hash, num, k, l);
	} else {
		*k = 1; *l = fmi->length - 1;
	}

	if (*k == 0 && *l == 0) return 0;

	return fmi_aln(fmi, qry, left, k, l);
}

void lc_write(const char *path, const lc_hash *h) {
	FILE *of = fopen(path, "w");
	fwrite(&h->hlen, sizeof(int), 1, of);
	fwrite(&h->len, sizeof(u64), 1, of);
	fwrite(h->lc, sizeof(u64), h->len, of);
	fclose(of);
}

void lc_read(const char *path, lc_hash *h) {
	FILE *inf = fopen(path, "r");
	int hlen;
	u64 len;
	u64 *lc;
	fread(&hlen, sizeof(int), 1, inf);
	fread(&len, sizeof(u64), 1, inf);
	lc = malloc(len * sizeof(u64));
	fread(lc, sizeof(u64), len, inf);
	h->len = len;
	h->hlen = hlen;
	h->lc = lc;
    fclose(inf);
}
