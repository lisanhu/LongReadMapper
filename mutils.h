//
// Created by lisanhu on 9/16/18.
//

#ifndef ACCSEQV8_MUTILS_H
#define ACCSEQV8_MUTILS_H

#include <stdint.h>
#include <time.h>
#include <stdio.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * mstring class
 */
typedef struct mstring {
    uint64_t l;
    char *s;
    bool own;
} mstring;

char * cstr_concat(const char *s1, const char *s2);
#pragma acc routine
mstring mstring_from(char *s, bool own);
#pragma acc routine
mstring mstring_borrow(char *s, size_t l);

#pragma acc routine
mstring mstring_own(const char *s, size_t l);
void mstring_destroy(mstring *ms);
void mstring_write(mstring ms, FILE *fp);
size_t mstring_read(mstring *ms, FILE *fp);


double time_elapse(struct timespec start);


size_t file_length(const char *path);

const char * load_file(const char *path, uint64_t *len);

char * cigar_align(const char *qry, int qlen, const char *target, int tlen, int *limit);

#ifdef __cplusplus
};
#endif

#endif //ACCSEQV8_MUTILS_H
