//
// Created by lisanhu on 9/16/18.
//

#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#include "mutils.h"
#include "edlib/edlib.h"
#include "gact/gact.h"

#pragma acc routine 
size_t acc_strlen(char * s) {
    size_t res = 0;
    while (*(s++) != '\0') res++;
    return res;
}

char * cstr_concat(const char *s1, const char *s2) {
    size_t l = strlen(s1) + strlen(s2);
    char *s = malloc(sizeof(char) * (l + 1));
    sprintf(s, "%s%s", s1, s2);
    s[l] = '\0';
    return s;
}

mstring mstring_from(char *s, bool own) {
    size_t l = acc_strlen(s);
    if (own) return mstring_own(s, l);
    else return mstring_borrow(s, l);
}


void mstring_destroy(mstring *ms) {
    if (ms->s && ms->own) {
        free(ms->s);
        ms->s = NULL;
    }
    ms->l = 0;
}


double time_elapse(struct timespec start) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return now.tv_sec - start.tv_sec +
           (now.tv_nsec - start.tv_nsec) / 1000000000.;
}


void mstring_write(mstring ms, FILE *fp) {
    fwrite(&ms.l, sizeof(uint64_t), 1, fp);
    fwrite(ms.s, sizeof(char), ms.l, fp);
}


size_t mstring_read(mstring *ms, FILE *fp) {
    size_t l = fread(&ms->l, sizeof(uint64_t), 1, fp);
    if (l == 0)
        return l;
    ms->s = malloc((ms->l + 1) * sizeof(char));
    l = fread(ms->s, sizeof(char), ms->l, fp);
    ms->s[ms->l] = '\0';
    ms->own = true;
	return l;
}

size_t file_length(const char *path) {
	FILE *fp = fopen(path, "r");
	fseek(fp, 0, SEEK_END);
	long l = ftell(fp);
	fclose(fp);

	if (l == EOF) {
		// todo: error case
		return 0;
	} else {
		return (size_t)l;
	}
}

const char * load_file(const char *path, uint64_t *len) {
	FILE *fp = fopen(path, "r");
	*len = file_length(path);
	char *buf = malloc(*len + 1);
	buf[*len] = '\0';
	fread(buf, sizeof(char), *len, fp);
    fclose(fp);
	return buf;
}

inline cigar cigar_align(const char *qry, int qlen, const char *target,
                         int tlen, int *limit, uint8_t *cigar_result) { // NOLINT(readability-non-const-parameter)

    mmstring q = ms_borrow((char *) qry, qlen);
    mmstring d = ms_borrow((char *) target, tlen);
    cigar result ;//= {.cigar = cigar_result, .n_cigar_op = 0};
    result.cigar = cigar_result; 
    result.n_cigar_op = 0;
    simple_gact(q, d, &result);
    *limit = result.score;
    return result;
}

#pragma acc routine
mstring mstring_borrow(char *s, size_t l) { // NOLINT(readability-non-const-parameter)
    /// todo: report this problem to CLion or Clang
    /// why is this s could be const? ms contains the reference to s and
    /// can be used to modify s in the future
    mstring ms ;//= {.s = s, .l = l, .own = false};
    ms.s = s; 
    ms.l = l;
    ms.own = false;
    return ms;
}

mstring mstring_own(const char *s, size_t l) {
    mstring ms = {.s = strdup(s), .l = l, .own = true};
    return ms;
}

mstring mstring_clone(mstring ms) {
    mstring result = {.s = strndup(ms.s, ms.l), .l = ms.l, .own = true};
    return result;
}


