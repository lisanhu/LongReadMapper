cmake_minimum_required(VERSION 3.5)
project(AccSeqV8 C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(ZLIB REQUIRED)

add_subdirectory(psascan)
#add_subdirectory(edlib)
add_subdirectory(mlog)
add_subdirectory(gact)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
	# using GCC
    set(CMAKE_CXX_FLAGS "-fopenmp ${CMAKE_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "-fopenmp ${CMAKE_C_FLAGS}")

	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#		set(CMAKE_CXX_FLAGS "-pg -fsanitize=address ${CMAKE_CXX_FLAGS}")
		set(CMAKE_CXX_FLAGS "-pg ${CMAKE_CXX_FLAGS}")
#		set(CMAKE_EXE_LINKER_FLAGS "-pg -fsanitize=address ${CMAKE_EXE_LINKER_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "-pg ${CMAKE_EXE_LINKER_FLAGS}")
	elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
		set(CMAKE_C_FLAGS "-O3 ${CMAKE_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "-O3 ${CMAKE_CXX_FLAGS}")
	endif ()
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
	# using PGI
	# todo: PGI bug work-around
#	set(CMAKE_CXX_FLAGS_DEBUG "-O0")
#	set(CMAKE_CXX_FLAGS_RELEASE "-O0 -g ")
#	set(CMAKE_C_FLAGS_RELEASE  "-O0 -g ")
#	set(CMAKE_C_FLAGS_DEBUG "-pg ${CMAKE_C_FLAGS_DEBUG}")
	#    message("${CMAKE_CXX_FLAGS_DEBUG}")

	#set(ACCTYPE "multicore")

	if (NOT ACCTYPE)
		set(ACCTYPE "s")
	endif ()
	
	# Adding ACC flags 
#	if (NOT ("${ACCTYPE}" STREQUAL "s" OR "${ACCTYPE}" STREQUAL "serial"))
#		if (OMPMULTICORE AND "${ACCTYPE}" STREQUAL "multicore")
#			message (WARNING "ACC Multicore and OMP multicore conflict. ACC Has priority")
#		endif ()
#		set(CMAKE_CXX_FLAGS "-ta=${ACCTYPE} -Minfo ${CMAKE_CXX_FLAGS}")
#		set(CMAKE_C_FLAGS "-ta=${ACCTYPE} -Minfo ${CMAKE_C_FLAGS}")
#	endif ()
	# Adding OMP Flags	
	if (OMPMULTICORE AND NOT "${ACCTYPE}" STREQUAL "multicore")
		set(CMAKE_CXX_FLAGS "-mp ${CMAKE_CXX_FLAGS}")
		set(CMAKE_C_FLAGS "-mp ${CMAKE_C_FLAGS}")
	endif ()

	# Selecting multicore flags
	if ("${ACCTYPE}" STREQUAL "multicore")
		set(CMAKE_CXX_FLAGS "-DMP_PARALLELISM=1 ${CMAKE_CXX_FLAGS}")
	elseif (OMPMULTICORE)
		set(CMAKE_CXX_FLAGS "-DMP_PARALLELISM=2 ${CMAKE_CXX_FLAGS}")
	endif ()
	
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
	# using Clang
	set(CMAKE_CXX_FLAGS "-DMP_PARALLELISM=1 -fopenmp=libomp -g ${CMAKE_CXX_FLAGS}")
	set(CMAKE_C_FLAGS "-DMP_PARALLELISM=1 -fopenmp=libomp -g ${CMAKE_C_FLAGS}")

	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#		set(CMAKE_CXX_FLAGS "-pg -fsanitize=address ${CMAKE_CXX_FLAGS}")
		set(CMAKE_CXX_FLAGS "-pg ${CMAKE_CXX_FLAGS}")
#		set(CMAKE_EXE_LINKER_FLAGS "-pg -fsanitize=address ${CMAKE_EXE_LINKER_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "-pg ${CMAKE_EXE_LINKER_FLAGS}")
	elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
		set(CMAKE_C_FLAGS "-O3 ${CMAKE_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "-O3 ${CMAKE_CXX_FLAGS}")
	endif ()
endif()

add_executable(accidx asindex.c mutils.c mutils.h fmidx/fmidx.c
		fmidx/fmidx.h lchash/lchash.c lchash/lchash.h)
target_link_libraries(accidx sa_use divsufsort64 divsufsort gact)


add_executable(accaln accaln.c alnmain.c histo/histo.c mutils.c fmidx/fmidx.c lchash/lchash.c mlog/logger.c gact/gact.c gact/mstring/mstring.c)
target_link_libraries(accaln sa_use divsufsort64 divsufsort)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
	target_link_libraries(accidx -mp)
	target_link_libraries(accaln -mp)
endif()

if (ZLIB_FOUND)
	include_directories(${ZLIB_INCLUDE_DIRS})
	target_link_libraries(accidx ${ZLIB_LIBRARIES})
	target_link_libraries(accaln ${ZLIB_LIBRARIES})
endif (ZLIB_FOUND)


#add_executable(test test/play.c logger/logger.c ../mutils.c)
#target_link_libraries(test edlib)
