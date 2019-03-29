//
// Created by Sanhu Li on 06/14/18.
//

#ifndef ACCSEQ_LOGGER_H
#define ACCSEQ_LOGGER_H

#include <stdio.h>

typedef struct timespec ts_t;
typedef struct mlog mlog;

typedef ts_t (*ml_mp_ptr)(mlog *, int, const char *, ...);
typedef ts_t (*ml_mv_ptr)(mlog *, const char *, ...);
typedef ts_t (*ml_me_ptr)(mlog *, const char *, ...);

struct mlog {
    ts_t start;
    FILE *stream;
    ml_mp_ptr mprint;
    ml_mv_ptr mvlog;
    ml_me_ptr melog;
};


mlog new_mlogger(ts_t *start);

#endif //ACCSEQ_LOGGER_H
