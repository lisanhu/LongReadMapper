//
// Created by Sanhu Li on 08/18/18.
//

#ifndef ACCSEQ_SA_USE_H
#define ACCSEQ_SA_USE_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>


#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	uint32_t    low;
	uint8_t     high;
} ui40_t;

static inline uint64_t ui40_convert(ui40_t val) {
	return ((uint64_t)val.high) << 32 | (uint64_t)val.low;
}

static inline size_t ui40_fread(ui40_t *buf, size_t nitems, FILE *stream) {
	size_t i = 0;
	uint32_t low;
	uint8_t high;
	printf("%s:%d\n", __FILE__, __LINE__);
	for (; i < nitems; ++i) {
		printf("%s:%d\n", __FILE__, __LINE__);
		fread(&low, sizeof(low), 1, stream);
		printf("%s:%d\n", __FILE__, __LINE__);
		size_t sz = fread(&high, sizeof(high), 1, stream);
		printf("%s:%d\n", __FILE__, __LINE__);
		if (sz == 0) break;
		buf[i].low = low;
		buf[i].high = high;
	}
	printf("%s:%d\n", __FILE__, __LINE__);
	return i;
}

void sa_build(const char *fname, long ram_use);

#ifdef __cplusplus
}
#endif


#endif //ACCSEQ_SA_USE_H
