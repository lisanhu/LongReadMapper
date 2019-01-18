//
// Created by lisanhu on 11/9/18.
//

#include <stdint.h>

#ifndef ACCSEQV8_SSW_USE_H
#define ACCSEQV8_SSW_USE_H


#ifdef __cplusplus
extern "C" {
#endif

char *compute_cigar(const char *s1, uint32_t l1, const char *s2, uint32_t l2, uint16_t *score);

#ifdef __cplusplus
};
#endif


#endif //ACCSEQV8_SSW_USE_H
