//
// Created by lisanhu on 9/17/18.
//
#include <stdio.h>
#include <time.h>
#include <zconf.h>
#include "../logger/logger.h"

int main(int argc, const char **argv) {
    unsigned int delay = 2;
    clock_gettime(CLOCK_MONOTONIC, &start);
    mprint(AS_LOG_VERBOSE, "I want to %s something.", "say");
    sleep(delay);
    mprint(AS_LOG_ERROR, "This is an %s message that delayed %d seconds...", "Error", delay);
	return 0;
}