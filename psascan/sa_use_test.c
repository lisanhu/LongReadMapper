//
// Created by Sanhu Li on 08/18/18.
//

#include <stdio.h>
#include "sa_use.h"

int main() {
	const char *fname = "input2.txt.sa5";
	const int BUF_LEN = 1000;
	ui40_t buf[BUF_LEN];


	sa_build("input2.txt", 8L << 30);


	FILE *in = fopen(fname, "r");
	size_t len = ui40_fread(buf, BUF_LEN, in);
	for (int i = 0; i < len; ++i) {
		printf("%lu\n", ui40_convert(buf[i]));
	}

	return 0;
}


