//
// Created by lisanhu on 3/29/19.
//

#include <stdlib.h>
#include <unistd.h>
#include "logger.h"

int main(int argc, const char *argv[]) {
    mlog logger = new_mlogger(NULL);
    logger.mvlog(&logger, "This is a verbose log without any dynamic values.");
    logger.melog(&logger, "This is an %s log with dynamic values.", "Error");
    logger.mvlog(&logger, "%s %s %s %s %s %s %s %s %s %s!", "You", "can",
                 "actually", "use", "printf", "format", "and", "argument",
                 "rules", "here");
    unsigned int delay = 1;
    sleep(delay);
    logger.mvlog(&logger, "You can definitely delay for %d second(s) to see "
                          "the timer on the left!", delay);
    logger.stream = stdout;
    logger.mvlog(&logger,
                 "By default, we are writing log info to stderr to make it "
                 "easier to separate your log information and output. But you "
                 "can definitely change the output stream any time you want!");
    return 0;
}