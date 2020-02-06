cmake_minimum_required(VERSION 3.5)
project(AccSeqV8 C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE ON)

# find_package(ZLIB REQUIRED)

add_subdirectory(psascan)
add_subdirectory(edlib)
add_subdirectory(mlog)

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
	set(CMAKE_CXX_FLAGS_DEBUG "-O0")
	set(CMAKE_CXX_FLAGS_RELEASE "-O4 ${CMAKE_CXX_FLAGS_RELEASE}")
	set(CMAKE_C_FLAGS_DEBUG "-pg ${CMAKE_C_FLAGS_DEBUG}")
	set(CMAKE_C_FLAGS_RELEASE  "-O4 ${CMAKE_C_FLAGS_RELEASE}")
	#    message("${CMAKE_CXX_FLAGS_DEBUG}")

	set(ACCTYPE "multicore")

	if (NOT ACCTYPE)
		set(ACCTYPE "s")
	endif ()

	if (NOT ("${ACCTYPE}" STREQUAL "s" OR "${ACCTYPE}" STREQUAL "serial"))
		set(CMAKE_CXX_FLAGS "-ta=${ACCTYPE} -Minfo=accel ${CMAKE_CXX_FLAGS}")
		set(CMAKE_C_FLAGS "-ta=${ACCTYPE} -Minfo=accel ${CMAKE_C_FLAGS}")
	endif ()
	
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
	# using Clang
	set(CMAKE_CXX_FLAGS "-fopenmp=libomp -g ${CMAKE_CXX_FLAGS}")
	set(CMAKE_C_FLAGS "-fopenmp=libomp -g ${CMAKE_C_FLAGS}")

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
target_link_libraries(accidx edlib sa_use divsufsort64 divsufsort z)


add_executable(accaln accaln.c alnmain.c histo/histo.c mutils.c fmidx/fmidx.c lchash/lchash.c)
target_link_libraries(accaln mlog edlib sa_use divsufsort64 divsufsort z)

# if (ZLIB_FOUND)
# 	include_directories(${ZLIB_INCLUDE_DIRS})
# 	target_link_libraries(accidx ${ZLIB_LIBRARIES})
# 	target_link_libraries(accaln ${ZLIB_LIBRARIES})
# endif (ZLIB_FOUND)


#add_executable(test test/play.c logger/logger.c ../mutils.c)
#target_link_libraries(test edlib)
