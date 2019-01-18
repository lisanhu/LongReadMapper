#include "sa_use.h"
#include "src/psascan_src/psascan.h"

#include <iostream>
#include <string>

void sa_build(const char *fname, long ram_use) {
	std::string t_fname(fname);
	std::string out_fname = t_fname + ".sa5";
	std::string gap_fname = out_fname;

//	long ram_use = 3072L << 20;
//	long ram_use = 8L << 30;
	long max_threads = (long)omp_get_max_threads();

	pSAscan(t_fname, out_fname, gap_fname,
	        ram_use, max_threads, false);
}
