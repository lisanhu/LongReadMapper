cmake_minimum_required(VERSION 3.5)
project(AccSeqV8 C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE ON)


add_subdirectory(psascan)
add_subdirectory(edlib)
#add_subdirectory(ssw)

if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
	set(CMAKE_CXX_FLAGS "-pg -fsanitize=address ${CMAKE_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "-pg -fsanitize=address ${CMAKE_EXE_LINKER_FLAGS}")
	message('test')
endif ()

#set(CMAKE_CXX_FLAGS "-pg ${CMAKE_CXX_FLAGS}")
#set(CMAKE_EXE_LINKER_FLAGS "-pg ${CMAKE_EXE_LINKER_FLAGS}")

add_executable(accidx asindex.c mutils.c mutils.h fmidx/fmidx.c fmidx/fmidx.h lchash/lchash.c lchash/lchash.h)
target_link_libraries(accidx -fopenmp)

add_executable(accaln accaln.c alnmain.c logger/logger.c histo/histo.c mutils.c fmidx/fmidx.c lchash/lchash.c)
target_link_libraries(accaln -fopenmp)
find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
	include_directories(${ZLIB_INCLUDE_DIRS})
	target_link_libraries(accidx ${ZLIB_LIBRARIES})
	target_link_libraries(accaln ${ZLIB_LIBRARIES})
endif (ZLIB_FOUND)
target_link_libraries(accidx sa_use edlib)
target_link_libraries(accaln sa_use edlib)


#add_executable(test test/test-main.cc test/test-fmidx.cc fmidx/fmidx.c mutils.c test/test-lchash.cc lchash/lchash.c)
#target_link_libraries(test sa_use -fopenmp)


#add_executable(play test/play.c mutils.c lchash/lchash.c fmidx/fmidx.c)
#target_link_libraries(play sa_use -fopenmp)