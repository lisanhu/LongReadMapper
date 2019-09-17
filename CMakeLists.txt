cmake_minimum_required(VERSION 3.5)
project(AccSeqV8 C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(ZLIB REQUIRED)

add_subdirectory(psascan)
add_subdirectory(mlog)
add_subdirectory(gact)
add_subdirectory(chaining)

add_executable(accidx asindex.c mutils.c mutils.h fmidx/fmidx.c
        fmidx/fmidx.h lchash/lchash.c lchash/lchash.h)
target_link_libraries(accidx z)
target_link_libraries(accidx sa_use divsufsort64 divsufsort gact)

add_executable(accaln accaln.c alnmain.c histo/histo.c mutils.c fmidx/fmidx.c
        lchash/lchash.c mlog/logger.c gact/gact.c gact/mstring/mstring.c)
target_link_libraries(accaln z sa_use divsufsort64 divsufsort chain)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    target_link_options(accidx PRIVATE -fopenmp)
    target_link_options(accaln PRIVATE -fopenmp)

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(accaln PRIVATE -pg -fsanitize=address)
        target_link_libraries(accaln asan)
        #	elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    endif ()
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    # using PGI
    # todo: PGI bug work-around
    target_compile_options(accaln PRIVATE -O0)
    target_link_options(accaln PRIVATE -mp)
    #target_link_options(accidx PRIVATE -mp)
    #	set(CMAKE_CXX_FLAGS_DEBUG "-O0")
    #	set(CMAKE_CXX_FLAGS_RELEASE "-O0 -g ")
    #	set(CMAKE_C_FLAGS_RELEASE  "-O0 -g ")
    #	set(CMAKE_C_FLAGS_DEBUG "-pg ${CMAKE_C_FLAGS_DEBUG}")
    #    message("${CMAKE_CXX_FLAGS_DEBUG}")

    #set(ACCTYPE "multicore")

    if (NOT ACCTYPE)
        set(ACCTYPE "serial")
    endif ()

    if ("${ACCTYPE}" STREQUAL "serial" OR "${ACCTYPE}" STREQUAL "s")
    elseif ("${ACCTYPE}" STREQUAL "mc" OR "${ACCTYPE}" STREQUAL "multicore")
        set(ACCTYPE "multicore")
        target_compile_options(accidx PRIVATE -ta=${ACCTYPE})
        target_compile_options(accaln PRIVATE -ta=${ACCTYPE})
    elseif ("${ACCTYPE}" STREQUAL "gpu" or "${ACCTYPE}" STREQUAL "nvidia")
        set(ACCTYPE "nvidia")
        target_compile_options(accidx PRIVATE -ta=${ACCTYPE})
        target_compile_options(accaln PRIVATE -ta=${ACCTYPE})
    endif ()
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    # using Clang
    target_compile_options(accidx PRIVATE -DMP_PARALLELISM=1 -g -fopenmp=libomp)
    target_compile_options(accaln PRIVATE -DMP_PARALLELISM=1 -g -fopenmp=libomp)
    target_link_options(accidx PRIVATE -fopenmp)
    target_link_options(accaln PRIVATE -fopenmp)

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(accaln PRIVATE -pg)
        target_compile_options(accidx PRIVATE -pg)
        #	elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    endif ()
endif ()
