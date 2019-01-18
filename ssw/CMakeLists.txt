cmake_minimum_required(VERSION 3.5)
project(SSW C)

add_library(ssw ssw.c ssw_use.c)

add_executable(ssw_test ssw_test.c)
target_link_libraries(ssw_test ssw)