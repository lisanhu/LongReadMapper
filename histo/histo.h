//
// Created by lisanhu on 9/30/18.
//

#ifndef ACCSEQV8_HISTO_H
#define ACCSEQV8_HISTO_H


#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdlib.h>

#include "khash.h"


typedef uint64_t u64;
typedef uint32_t u32;


typedef struct {
	u64 key, val;
} entry;

KHASH_MAP_INIT_INT64(64, u64)

typedef struct {
	entry *entries;
	u32 cap, size, maxi, smaxi; // max id and second max id
//	khiter_t k;
	khash_t(64) *kh;
} histo;


histo *histo_init(u32 cap);
void histo_add(histo *h, u64 key);
void histo_destroy(histo *h);
u64 histo_get(histo *h, u64 key);
u64 histo_find_max(histo *h);
u64 histo_find_2_max(histo *h, entry *store);

#ifdef __cplusplus
};
#endif


#endif //ACCSEQV8_HISTO_H
