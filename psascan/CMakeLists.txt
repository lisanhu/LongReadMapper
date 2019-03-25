cmake_minimum_required(VERSION 3.9)
project(pSAscan C CXX)

set(CMAKE_C_STANDARD 11)

add_subdirectory(libdivsufsort)


set(SRC_FILES src/psascan_src/utils.cpp src/main.cpp)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    # using Clang
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    # using PGI

    # todo: PGI bug work-around
    set(CMAKE_CXX_FLAGS_DEBUG "-O0")
#    message("${CMAKE_CXX_FLAGS_DEBUG}")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mp")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_FLAGS "-O3 ${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "-O3 ${CMAKE_CXX_FLAGS}")
endif ()

add_executable(psascan ${SRC_FILES})
target_include_directories(psascan PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/libdivsufsort/include)
target_link_libraries(psascan divsufsort divsufsort64)


add_library(sa_use sa_use.cc src/psascan_src/utils.cpp)
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    target_compile_options(sa_use BEFORE PRIVATE "-fopenmp")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    # using PGI
    target_compile_options(sa_use BEFORE PRIVATE "-mp")
endif()

target_include_directories(sa_use PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/libdivsufsort/include)
target_link_libraries(sa_use divsufsort divsufsort64)

add_executable(sa_use_test sa_use_test.c)
target_link_libraries(sa_use_test sa_use)
