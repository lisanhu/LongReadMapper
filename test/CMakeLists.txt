cmake_minimum_required(VERSION 3.5)
project(tests C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE ON)

#add_subdirectory(../)


include_directories(../)
add_executable(accseq-tests test-main.cc test-fmidx.cc ../fmidx/fmidx.c ../mutils.c ../gact/mutils.c)
target_link_libraries(accseq-tests gact sa_use)
