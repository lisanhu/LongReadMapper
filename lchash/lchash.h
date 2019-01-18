//
// Created by lisanhu on 9/24/18.
//

#ifndef ACCSEQV8_LCHASH_H
#define ACCSEQV8_LCHASH_H

#include <stdint.h>

#include "../fmidx/fmidx.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	uint64_t *lc;
	uint64_t len;
	int hlen;
} lc_hash;

uint64_t lc_access(lc_hash h, uint64_t loc, uint64_t *k, uint64_t *l);
lc_hash lc_build(const dna_fmi *fmi, int hlen);
uint64_t lc_length(int hlen);
uint64_t lc_aln(const char *qry, int qlen, uint64_t *k, uint64_t *l,
                const dna_fmi *fmi,
                const lc_hash *hash);
void lc_destroy(lc_hash hash);
void lc_write(const char *path, const lc_hash *h);
void lc_read(const char *path, lc_hash *h);

#ifdef __cplusplus
};
#endif


#endif //ACCSEQV8_LCHASH_H
