//
// Created by lisanhu on 11/9/18.
//

#include <stdio.h>
#include <string.h>
#include "ssw_use.h"

int main() {
	uint16_t s;
	const char *a = "GAATTC";
	const char *b = "GAATTCC";
	uint32_t l1 = (uint32_t) strlen(a);
	uint32_t l2 = (uint32_t) strlen(b);

	compute_cigar(a, l1, b, l2, &s);
	return 0;
}