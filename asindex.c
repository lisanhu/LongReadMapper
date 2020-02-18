#include <stdio.h>
#include <zlib.h>

#include "fmidx/fmidx.h"
#include "kseq.h"
#include "lchash/lchash.h"
#include "mutils.h"
#include "psascan/sa_use.h"

KSEQ_INIT(gzFile, gzread)

char revc_mapper[256];

static inline void _seq_to_upper_case(char *seq, size_t length);

static mstring ms_from_ks(kstring_t ks) {
    mstring ms = {.l = ks.l, .s = strndup(ks.s, ks.l)};
    return ms;
}

static void _cstr_reverse_inplace(char *s) {
    size_t l = strlen(s);
    char c;
    for (size_t i = 0; i < l / 2; ++i) {
        c = s[i];
        s[i] = s[l - i - 1];
        s[l - i - 1] = c;
    }
}

static int32_t _dna_rand_ch() {
    static int32_t val = 0;
    static int pos = -1;
    if (pos < 0) {
        val = (int32_t)lrand48();
        pos = 0;
    } else if (pos < 31) {
        pos += 2;
    } else {
        val = (int32_t)lrand48();
        pos = 0;
    }
    return (val >> pos) & 0x3;
}

static bool charin(const char *seq, char c) {
    size_t l = strlen(seq);
    for (size_t i = 0; i < l; i++) {
        if (c == seq[i]) {
            return true;
        }
    }
    return false;
}

static void _dna_replace_n_inplace(char *seq, size_t length) {
    for (size_t i = 0; i < length; ++i) {
        if (!charin("GACTgact", seq[i])) {
            int c = _dna_rand_ch();
            seq[i] = "ACGT"[c];
        }
    }
}

static void _seq_to_upper_case(char *seq, size_t length) {
    for (size_t i = 0; i < length; ++i) {
        if (seq[i] > 0x60)
            seq[i] -= 0x20;
    }
}

static void _dna_rev_complementary_inplace(char *seq, size_t length) {
    for (size_t i = 0; i < length; ++i) {
        seq[i] = revc_mapper[seq[i]];
    }
    _cstr_reverse_inplace(seq);
}

void create_meta(const char *prefix) {
    char *mfn = cstr_concat(prefix, ".mta");
    FILE *mfp = fopen(mfn, "w");
    char *cfn = cstr_concat(prefix, ".cat");
    FILE *cfp = fopen(cfn, "w");

    uint64_t offset = 0;

    gzFile fp = gzopen(prefix, "r");
    kseq_t *seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        mstring seq_name = ms_from_ks(seq->name);
        /// mta file contents
        mstring_write(seq_name, mfp);                    // seq_name
        fwrite(&offset, sizeof(uint64_t), 1, mfp);       // offset
        fwrite(&seq->seq.l, sizeof(seq->seq.l), 1, mfp); // seq_len
        mstring_destroy(&seq_name);

        /// cat file contents
        char *seqs = strndup(seq->seq.s, seq->seq.l);
        _dna_replace_n_inplace(seqs, seq->seq.l);
        _seq_to_upper_case(seqs, seq->seq.l);
        offset += fwrite(seqs, sizeof(char), seq->seq.l, cfp);
        _dna_rev_complementary_inplace(seqs, seq->seq.l);
        offset += fwrite(seqs, sizeof(char), seq->seq.l, cfp);
        free(seqs);
    }
    kseq_destroy(seq);
    gzclose(fp);

    /// put a '$' at the end of cat file
    char c = '$';
    fwrite(&c, sizeof(char), 1, cfp);

    fclose(cfp);
    free(cfn);
    fclose(mfp);
    free(mfn);
}

void init() {
    revc_mapper['a'] = revc_mapper['A'] = 'T';
    revc_mapper['c'] = revc_mapper['C'] = 'G';
    revc_mapper['g'] = revc_mapper['G'] = 'C';
    revc_mapper['t'] = revc_mapper['T'] = 'A';

    srand48(time(NULL));
}

int main(int argc, const char **argv) {
    init();

    create_meta(argv[1]);
    char *cfn = cstr_concat(argv[1], ".cat");

    dna_fmi idx;
    //	fmi_build(cfn, &idx, 128, 8L << 30);
    printf("Start building fmindex\n");
    fmi_build(cfn, &idx, 32, 8L << 30);
    printf("Done building fmindex, start writing to disk\n");
    fmi_write(idx, cfn);

    int hash_len = 12;
    printf("Start building lchash\n");
    lc_hash hash = lc_build(&idx, hash_len);
    char *lch_path = cstr_concat(cfn, ".lch");
    printf("Start writing lchash to disk\n");
    lc_write(lch_path, &hash);
    lc_destroy(hash);

    free(lch_path);
    free(cfn);
    return 0;
}
