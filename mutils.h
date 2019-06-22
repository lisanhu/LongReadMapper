//
// Created by lisanhu on 9/16/18.
//

#ifndef ACCSEQV8_MUTILS_H
#define ACCSEQV8_MUTILS_H

#define ACC_PARALLEL 1
#define OMP_PARALLEL 2

#ifndef MP_PARALLELISM
#define MP_PARALLELISM OMP_PARALLEL
#endif


#if MP_PARALLELISM == OMP_PARALLEL
#include <omp.h>
#endif
#include <stdint.h>
#include <time.h>
#include <stdio.h>
#include <stdbool.h>
#include "gact/gact.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * mstring class
 */
typedef struct mstring {
    uint64_t l;
    char *s;
    int own;
} mstring;

char * cstr_concat(const char *s1, const char *s2);
mstring mstring_from(char *s, bool own);
#pragma acc routine
mstring mstring_borrow(char *s, size_t l);
mstring mstring_own(const char *s, size_t l);
mstring mstring_clone(mstring ms);
void mstring_destroy(mstring *ms);
void mstring_write(mstring ms, FILE *fp);
size_t mstring_read(mstring *ms, FILE *fp);


double time_elapse(struct timespec start);


size_t file_length(const char *path);

const char * load_file(const char *path, uint64_t *len);

#pragma acc routine 
cigar cigar_align(const char *qry, int qlen, const char *target, int tlen,
                  int *limit, uint8_t *cigar_result);

#ifdef __cplusplus
};
#endif

#endif //ACCSEQV8_MUTILS_H
