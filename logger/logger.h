//
// Created by Sanhu Li on 06/14/18.
//

#ifndef ACCSEQ_LOGGER_H
#define ACCSEQ_LOGGER_H

#define AS_LOG_VERBOSE  0
#define AS_LOG_ERROR    1

typedef struct timespec ts_t;

extern int level;

struct timespec print_log(int l, const char *text, ts_t start_time);

#endif //ACCSEQ_LOGGER_H
