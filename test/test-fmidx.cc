//
// Created by lisanhu on 9/21/18.
//

#include <cstring>
#include <cstdlib>

#include "catch.hpp"
#include "fmidx/fmidx.h"
#include "psascan/sa_use.h"
#include "mutils.h"

TEST_CASE("fmindex tests", "[fmidx]") {
	dna_fmi fmi;
	fmi_build("input2.txt", &fmi, 128, 8L<<30);
	fmi_write(fmi, "input2.txt");

	SECTION("write and loading mfi file from disk") {
		dna_fmi mfi;
		fmi_read(&mfi, "input2.txt");
		REQUIRE(fmi.o_ratio == mfi.o_ratio);
		REQUIRE(fmi.length == mfi.length);
		REQUIRE(!memcmp(fmi.bwt, mfi.bwt, fmi.length));
		REQUIRE(fmi.o_len == mfi.o_len);
		REQUIRE(!memcmp(fmi.o, mfi.o, fmi.o_len * sizeof(uint64_t)));
		REQUIRE(!memcmp(fmi.c, mfi.c,256));
		fmi_destroy(&mfi);
	}

	SECTION("perform aln") {
		const char *query = "atgcgcatatgtgtaaaacgtaactagatgtcgggaaggcggagcaagagcagctacatcaatg";
		int len = static_cast<int>(strlen(query));
		uint64_t k = 1, l = fmi.length - 1;
		uint64_t r = fmi_aln(&fmi, query, len, &k, &l);
		REQUIRE(r == 1);
		uint64_t off = sa_access("input2.txt", 1024, k);
		uint64_t flen;
		const char *content = load_file("input2.txt", &flen);
		REQUIRE(!memcmp(query, content + off, len));
		delete[] content;
	}

	fmi_destroy(&fmi);
}


/**
 * The strange output is a number 1 appended after all the input, like:
 * ```
 * blablabla.... #correct output
 * 1 # in error stream
 * ```
 * Reason is catch2 output a xml report (used by CLion)
 * The output contains <StdErr>1</StdErr>
 * CLion seems to interpret it as a stderr output
 *
 * Probably this is because we have write operation for stderr, and is
 * captured by catch2
 *
 * Need to test with the 2018.3 EAP to see whether there is improved behaviour
 */
TEST_CASE("sa_build strange output", "[sa_build]") {
	FILE *err = stderr;
	FILE *nfp = fopen("/dev/null", "w");
	stderr = nfp;
	sa_build("input2.txt", 2L << 30);
	fclose(nfp);
	stderr = err;
}

static void  _build_csa(dna_fmi *idx, int *occ, int *c, const int *sa, size_t l);

TEST_CASE("Test CSA", "[sa_build]") {
    printf("Testing CSA\n");
    dna_fmi fmi;
    fmi.bwt = strdup("ard$rcaaaabb");
    int c_tab[128];
    size_t l = strlen(fmi.bwt);
    c_tab['$'] = 0;
    c_tab['a'] = 1;
    c_tab['b'] = 6;
    c_tab['c'] = 8;
    c_tab['d'] = 9;
    c_tab['r'] = 10;

    int o_tab[72] = {0, 1, 0, 0, 0, 0,
                     0, 1, 0, 0, 0, 1,
                     0, 1, 0, 0, 1, 1,
                     1, 1, 0, 0, 1, 1,
                     1, 1, 0, 0, 1, 2,
                     1, 1, 0, 1, 1, 2,
                     1, 2, 0, 1, 1, 2,
                     1, 3, 0, 1, 1, 2,
                     1, 4, 0, 1, 1, 2,
                     1, 5, 0, 1, 1, 2,
                     1, 5, 1, 1, 1, 2,
                     1, 5, 2, 1, 1, 2};

    int sa[] = {11, 10, 7, 0, 3, 5, 8, 1, 4, 6, 9, 2};

    printf("Start building CSA\n");
    fmi.csa_ratio = 32;
    _build_csa(&fmi, o_tab, c_tab, sa, l);
    for (int i = 0; i < l; ++i) {
        uint64_t si = csa_access(&fmi, i);
        CHECK(si == sa[i]);
    }
}

void _build_csa(dna_fmi *idx, int *occ, int *c, const int *sa, size_t l) {
    size_t csa_l = l / idx->csa_ratio + 1;
    idx->csa = static_cast<uint64_t *>(malloc(sizeof(uint64_t) * csa_l));

    for (uint64_t i = 0; i < csa_l; ++i) {
        uint64_t si = sa[i * idx->csa_ratio];
        idx->csa[i] = si;
    }
}
