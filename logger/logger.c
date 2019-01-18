//
// Created by Sanhu Li on 06/14/18.
//

#include <stdio.h>
#include <time.h>
#include "logger.h"
#include "../mutils.h"

int level = AS_LOG_VERBOSE;

struct timespec print_log(int l, const char *text, ts_t start_time) {
	if (l >= level) {
		char *ltext;
		switch (l) {
			case AS_LOG_VERBOSE:
				ltext = "Verbose";
				break;
			case AS_LOG_ERROR:
				ltext = "Error";
				break;
			default:
				ltext = "UNKNOWN";
				break;
		}

		ts_t tp;
		clock_gettime(CLOCK_MONOTONIC, &tp);
		double elapsed = time_elapse(start_time);

		fprintf(stderr, "[ %-7s %8.2lf ] %s\n", ltext, elapsed, text);
		return tp;
	}
}
