//
// Created by lisanhu on 3/20/19.
//

#ifndef ACCSEQV8_ALNMAIN_H
#define ACCSEQV8_ALNMAIN_H

#include "accaln.h"

typedef struct params {
	u64 batch_size;
	u32 seed_len, thres;
} params;

params read_params(const char *path, context *ctx);

#endif //ACCSEQV8_ALNMAIN_H
