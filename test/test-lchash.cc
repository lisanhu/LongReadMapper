//
// Created by lisanhu on 9/21/18.
//

#include <string.h>

#include "catch.hpp"
#include "../fmidx/fmidx.h"
#include "../mutils.h"
#include "../lchash/lchash.h"
#include "../accaln.h"

TEST_CASE("lchash tests", "[lchash]") {
	dna_fmi fmi;
	const char *prefix = "input2-cap.txt";
	fmi_build(prefix, &fmi, 128, 8L<<30);
	lc_hash h = lc_build(&fmi, 12);

	SECTION("lchash serialization") {
		char *path = cstr_concat(prefix, ".lch");
		lc_write(path, &h);
		lc_hash h2;
		lc_read(path, &h2);
		REQUIRE(h.hlen == h2.hlen);
		REQUIRE(h.len == h2.len);
		REQUIRE(!memcmp(h.lc, h2.lc, h.len * sizeof(u64)));
		free(path);
	}

	SECTION("lchash alignment") {
//		const char *query = "atgcgcatatgtgtaaaacgtaactagatgtcgggaaggcggagcaagagcagctacatcaatg";
		const char *query = "ATGCGCATATGTGTAAAACGTAACTAGATGTCGGGAAGGCGGAGCAAGAGCAGCTACATCAATG";
		const char *dbg = "GCTACATCAATG";
		int len = static_cast<int>(strlen(query));
		uint64_t k = 1, l = fmi.length - 1;
		uint64_t r = fmi_aln(&fmi, query, len, &k, &l);
		uint64_t kk, ll, rr;
		rr = lc_aln(query, len, &kk, &ll, &fmi, &h);
		REQUIRE(r == 1);
		REQUIRE(rr == r);
		REQUIRE(kk == k);
		REQUIRE(ll == l);
		uint64_t off = sa_access(prefix, 1024, kk);
		uint64_t flen;
		const char *content = load_file(prefix, &flen);
		REQUIRE(!memcmp(query, content + off, len));
		delete[] content;
	}


	fmi_destroy(&fmi);
	lc_destroy(h);
}
