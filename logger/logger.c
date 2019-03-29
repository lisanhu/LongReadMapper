//
// Created by Sanhu Li on 06/14/18.
//

#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include "logger.h"


#define AS_LOG_VERBOSE  0
#define AS_LOG_ERROR    1


static ts_t vmprint(mlog *tthis, int l, const char *fmt, va_list args);
static ts_t mvlog(mlog *tthis, const char *fmt, ...);
static ts_t melog(mlog *tthis, const char *fmt, ...);
static ts_t mprint(mlog *tthis, int l, const char *fmt, ...);
static inline double time_elapse(struct timespec start);


mlog new_mlogger(ts_t *start) {
    mlog logger;
    logger.stream = stderr;
    if (start == NULL) {
        clock_gettime(CLOCK_MONOTONIC, &logger.start);
    } else {
        logger.start = *start;
    }
    logger.mprint = mprint;
    logger.mvlog = mvlog;
    logger.melog = melog;
    return logger;
}


ts_t mvlog(mlog *tthis, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    ts_t tp = vmprint(tthis, AS_LOG_VERBOSE, fmt, args);
    va_end(args);
    return tp;
}

ts_t mprint(mlog *tthis, int l, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    ts_t tp = vmprint(tthis, l, fmt, args);
    va_end(args);
    return tp;
}

static ts_t vmprint(mlog *tthis, int l, const char *fmt, va_list args) {
    const char *ltext;
    const int LOG_BUF_SZ = 1024;
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
    char text[LOG_BUF_SZ];
    clock_gettime(CLOCK_MONOTONIC, &tp);
    double elapsed = time_elapse(tthis->start);

    vsprintf(text, fmt, args);

    fprintf(tthis->stream, "[ %-7s %8.2lf ] %s\n", ltext, elapsed, text);
    return tp;
}

ts_t melog(mlog *tthis, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    ts_t tp = vmprint(tthis, AS_LOG_ERROR, fmt, args);
    va_end(args);
    return tp;
}

static inline double time_elapse(struct timespec start) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return now.tv_sec - start.tv_sec +
           (now.tv_nsec - start.tv_nsec) / 1000000000.;
}