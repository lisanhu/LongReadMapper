//
// Created by lisanhu on 9/21/18.
//

#include <cstring>
#include <cstdlib>

#include "catch.hpp"
#include "../fmidx/fmidx.h"
#include "../psascan/sa_use.h"
#include "../mutils.h"

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
