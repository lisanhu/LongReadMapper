//
// Created by lisanhu on 9/16/18.
//

#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "mutils.h"
#include "edlib/edlib.h"


char * cstr_concat(const char *s1, const char *s2) {
    size_t l = strlen(s1) + strlen(s2);
    char *s = malloc(sizeof(char) * (l + 1));
    sprintf(s, "%s%s", s1, s2);
    s[l] = '\0';
    return s;
}


mstring mstring_from(char *s) {
    mstring ms = {.s = strdup(s), .l = strlen(s)};
    return ms;
}


void mstring_destroy(mstring *ms) {
    if (ms->s) {
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
    fread(&ms->l, sizeof(uint64_t), 1, fp);
    ms->s = malloc((ms->l + 1) * sizeof(char));
    size_t l = fread(ms->s, sizeof(char), ms->l, fp);
    ms->s[ms->l] = '\0';
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
	return buf;
}

inline char *cigar_align(const char *qry, int qlen, const char *target, int tlen,
                  int *limit) {
	EdlibAlignResult align = edlibAlign(qry, qlen,target, tlen,
			edlibNewAlignConfig(*limit, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
	char *cigar = edlibAlignmentToCigar(align.alignment, align.alignmentLength,
	                                    EDLIB_CIGAR_STANDARD);
	*limit = align.editDistance;
	edlibFreeAlignResult(align);
	return cigar;
}


