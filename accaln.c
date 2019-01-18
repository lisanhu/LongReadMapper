//
// Created by lisanhu on 9/26/18.
//
#include "accaln.h"

void context_destroy(context *ctx)  {
	if (ctx->read1 != NULL) {
		free(ctx->read1);
		ctx->read1 = NULL;
	}

	if (ctx->read2 != NULL) {
		free(ctx->read1);
		ctx->read1 = NULL;
	}
}

size_t reads_load(read_t *buf, size_t num, kseq_t *seq) {
	size_t i;
	for (i = 0; i < num && kseq_read(seq) >= 0; ++i) {
		u64 len = seq->seq.l;
		read_t r = {.rid = i + 1,
			  .name = mstring_from(seq->name.s),
			  .seq = strndup(seq->seq.s, len),
			  .len = (u32) len,
			  .qual = strndup(seq->qual.s, len)};
			  /// rid start with 1, so we can assume seed with rid=0 is end
		buf[i] = r;
	}
	return i;
}

void read_destroy(read_t *rd) {
	if (rd != NULL) {
		rd->rid = 0;
		mstring_destroy(&rd->name);
		MFREE(rd->seq)
		MFREE(rd->qual)
		rd->len = 0;
	}
}

void reads_destroy(read_t *rd, size_t num) {
	for (size_t i = 0; i < num; ++i) {
		mstring_destroy(&rd[i].name);
		free(rd[i].qual);
		free(rd[i].seq);
	}
	memset(rd, 0, num * sizeof(read_t));
}

