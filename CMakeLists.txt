cmake_minimum_required(VERSION 3.5)
project(AccSeqV8 C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(ZLIB REQUIRED)

add_subdirectory(psascan)
add_subdirectory(edlib)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
	# using GCC
    set(CMAKE_CXX_FLAGS "-fopenmp ${CMAKE_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "-fopenmp ${CMAKE_CXX_FLAGS}")

	if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
		set(CMAKE_CXX_FLAGS "-pg -fsanitize=address ${CMAKE_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "-pg fsanitize=address ${CMAKE_EXE_LINKER_FLAGS}")
	endif ()
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
	# using PGI
	# todo: PGI bug work-around
	set(CMAKE_CXX_FLAGS_DEBUG "-O0")
	#    message("${CMAKE_CXX_FLAGS_DEBUG}")

	set(CMAKE_CXX_FLAGS "-ta=multicore -Minfo ${CMAKE_CXX_FLAGS}")
	set(CMAKE_C_FLAGS "-ta=multicore -Minfo ${CMAKE_CXX_FLAGS}")

	if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
		set(CMAKE_CXX_FLAGS "-pg ${CMAKE_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "-pg ${CMAKE_EXE_LINKER_FLAGS}")
	endif ()
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
	# using Clang
	set(CMAKE_CXX_FLAGS "-fopenmp=libomp ${CMAKE_CXX_FLAGS}")
	set(CMAKE_C_FLAGS "-fopenmp=libomp ${CMAKE_CXX_FLAGS}")

	if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
		set(CMAKE_CXX_FLAGS "-pg -fsanitize=address ${CMAKE_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "-pg fsanitize=address ${CMAKE_EXE_LINKER_FLAGS}")
	endif ()
endif()

add_executable(accidx asindex.c mutils.c mutils.h fmidx/fmidx.c
		fmidx/fmidx.h lchash/lchash.c lchash/lchash.h)
target_link_libraries(accidx edlib sa_use divsufsort64 divsufsort)


add_executable(accaln accaln.c alnmain.c logger/logger.c
		histo/histo.c mutils.c fmidx/fmidx.c lchash/lchash.c)
target_link_libraries(accaln edlib sa_use divsufsort64 divsufsort)

if (ZLIB_FOUND)
	include_directories(${ZLIB_INCLUDE_DIRS})
	target_link_libraries(accidx ${ZLIB_LIBRARIES})
	target_link_libraries(accaln ${ZLIB_LIBRARIES})
endif (ZLIB_FOUND)