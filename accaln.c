//
// Created by lisanhu on 9/26/18.
//
#include "accaln.h"
#include "fmidx/fmidx.h"

void context_destroy(context *ctx)  {
#define TRY_FREE(x) if ((x) != NULL) {free((x)); (x) = NULL;}
//	if (ctx->read1 != NULL) {
//		free(ctx->read1);
//		ctx->read1 = NULL;
//	}
	TRY_FREE(ctx->read1)

//	if (ctx->read2 != NULL) {
//		free(ctx->read1);
//		ctx->read1 = NULL;
//	}
    TRY_FREE(ctx->read2)
    TRY_FREE(ctx->genome) // using strdup, strndup
    TRY_FREE(ctx->prefix) // using cstr_concat

	if (ctx->fmi != NULL) {
        fmi_destroy(ctx->fmi);
        TRY_FREE(ctx->fmi)
	}

	if (ctx->lch != NULL) {
        lc_destroy(*ctx->lch);
        TRY_FREE(ctx->lch)
	}

    for (int i = 0; i < ctx->mta_len; ++i) {
        mstring ms = ctx->mta[i].seq_name;
        mstring_destroy(&ms);   // using mstring_read
    }
	TRY_FREE(ctx->mta)
	TRY_FREE(ctx->content) // using load_file

    TRY_FREE(sa_buf->mem)
    TRY_FREE(sa_buf)
#undef TRY_FREE
}

size_t reads_load(read_t *buf, size_t num, kseq_t *seq) {
	size_t i;
	for (i = 0; i < num && kseq_read(seq) >= 0; ++i) {
		u64 len = seq->seq.l;
		read_t r = {.rid = i + 1,
			  .name = mstring_from(seq->name.s, true),
			  .seq = mstring_own(seq->seq.s, len),
			  .len = (u32) len,
			  .qual = mstring_own(seq->qual.s, len)};
			  /// rid start with 1, so we can assume seed with rid=0 is end
		buf[i] = r;
	}
	return i;
}

void read_destroy(read_t *rd) {
	if (rd != NULL) {
		rd->rid = 0;
		mstring_destroy(&rd->name);
        mstring_destroy(&rd->seq);
        mstring_destroy(&rd->qual);
		rd->len = 0;
	}
}

void reads_destroy(read_t *rd, size_t num) {
	for (size_t i = 0; i < num; ++i) {
		mstring_destroy(&rd[i].name);
		mstring_destroy(&rd[i].qual);
		mstring_destroy(&rd[i].seq);
	}
	memset(rd, 0, num * sizeof(read_t));
}

