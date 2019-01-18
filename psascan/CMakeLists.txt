cmake_minimum_required(VERSION 3.9)
project(pSAscan C CXX)

set(CMAKE_C_STANDARD 11)

add_subdirectory(libdivsufsort)


set(SRC_FILES src/psascan_src/utils.cpp src/main.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(psascan ${SRC_FILES})
target_include_directories(psascan PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/libdivsufsort/include)
target_link_libraries(psascan divsufsort divsufsort64)


add_library(sa_use sa_use.cc src/psascan_src/utils.cpp)
target_compile_options(sa_use BEFORE PRIVATE "-fopenmp")
target_include_directories(sa_use PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/libdivsufsort/include)
target_link_libraries(sa_use divsufsort divsufsort64)

add_executable(sa_use_test sa_use_test.c)
target_link_libraries(sa_use_test sa_use)
