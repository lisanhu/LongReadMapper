//
// Created by lisanhu on 11/9/18.
//

#include "ssw_use.h"
#include "ssw.h"
/* This table is used to transform nucleotide letters into numbers. */
static const int8_t nt_table[128] = {
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


// initialize scoring matrix for genome sequences
//  A  C  G  T	N (or other ambiguous code)
//  2 -2 -2 -2 	0	A
// -2  2 -2 -2 	0	C
// -2 -2  2 -2 	0	G
// -2 -2 -2  2 	0	T
//	0  0  0  0  0	N (or other ambiguous code)
static const int8_t s_mat[25] = {
		 2, -2, -2, -2, 0,
		-2,  2, -2, -2, 0,
		-2, -2,  2, -2, 0,
		-2, -2, -2,  2, 0,
		 0,  0,  0,  0, 0,
};


static inline int8_t *s_2_seq(const char *s, uint32_t l) {
	int8_t *seq = malloc(l * sizeof(int8_t));
	for (int i = 0; i < l; ++i) {
		seq[i] = nt_table[s[i]];
	}
	return seq;
}

char *
compute_cigar(const char *s1, uint32_t l1, const char *s2, uint32_t l2, uint16_t *score) {
	int8_t *seq1, *seq2;
//	const int32_t match = 2, mismatch = 2;
	const int32_t gap_open = 3, gap_extension = 1;

	seq1 = s_2_seq(s1, l1);
	seq2 = s_2_seq(s2, l2);

	s_profile *sp = ssw_init(seq1, l1, s_mat, 5, 2);
	s_align *sa = ssw_align (sp, seq2, l2, gap_open, gap_extension, 1, 0, 0, 15);

	align_destroy(sa);
	init_destroy(sp);
	free(seq2);
	free(seq1);

	printf("%d\n", sa->cigar[0]);

	for (int32_t i = 0; i < sa->cigarLen; ++i) {
		char letter = cigar_int_to_op(sa->cigar[i]);
		uint32_t length = cigar_int_to_len(sa->cigar[i]);
		printf("%d%c", length, letter);
	}
	printf("\n");

	return NULL;
}
