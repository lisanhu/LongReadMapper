//
// Created by lisanhu on 9/17/18.
//

#ifndef ACCSEQV8_FMIDX_H
#define ACCSEQV8_FMIDX_H


#include <stdint.h>
#include "../psascan/sa_use.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _dna_fmi {
	uint64_t length, o_len;
	uint64_t *c, *o;
	int o_ratio;
	char *bwt;
} dna_fmi;

typedef struct {
	uint64_t start, len;
	ui40_t *mem;
} sa_mem;

extern sa_mem *sa_buf;

uint64_t sa_access(const char *prefix, uint64_t cache_sz, uint64_t loc);
void sa_done_access();

void fmi_build(const char *prefix, dna_fmi *idx, int o_ratio, long ram_use);

uint64_t fmi_aln(const dna_fmi *idx, const char *qry, int len, uint64_t *k,
                 uint64_t *l);

void fmi_write(dna_fmi idx, const char *prefix);
void fmi_read(dna_fmi *idx, const char *prefix);

void fmi_destroy(dna_fmi *idx);


#ifdef __cplusplus
};
#endif

#endif //ACCSEQV8_FMIDX_H
