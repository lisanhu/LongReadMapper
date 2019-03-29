cmake_minimum_required(VERSION 3.5)

add_library(mlog logger.c)

option(BUILD_USE "Whether build the usage target" on)

if (BUILD_USE)
    add_executable(test_logger logger_use.c logger.c)
endif ()