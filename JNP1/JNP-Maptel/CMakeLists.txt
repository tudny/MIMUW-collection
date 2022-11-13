cmake_minimum_required(VERSION 3.0)
project(Zadanie2_maptel)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif ()

# We don't want to see make commands.
set(CMAKE_VERBOSE_MAKEFILE OFF)

# Settings flags for compiling C and CPP.

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-std=c++17")
set(CMAKE_C_FLAGS "-std=c11")

set(DEBUG_FLAGS "-g -Wall -Wextra")

set(CMAKE_CXX_FLAGS_DEBUG ${DEBUG_FLAGS})
set(CMAKE_C_FLAGS_DEBUG ${DEBUG_FLAGS})

set(RELEASE_FLAGS "-O2 -DNDEBUG")

set(CMAKE_CXX_FLAGS_RELEASE ${RELEASE_FLAGS})
set(CMAKE_C_FLAGS_RELEASE ${RELEASE_FLAGS})

# Settings source files for default executable file.

set(SOURCE_FILES
        src/test.c
        src/maptel.cc
        src/maptel.h )

add_executable(Zadanie2_maptel ${SOURCE_FILES})

# Setting tests

## Test for C compilation
set(SOURCE_TEST_1
        src/maptel.h
        src/maptel.cc
        src/maptel_test1.c)

add_executable(test1 EXCLUDE_FROM_ALL ${SOURCE_TEST_1})

set_target_properties(test1 PROPERTIES OUTPUT_NAME maptel_test1)

## Test for CPP compilation
set(SOURCE_TEST_2
        src/maptel.h
        src/maptel.cc
        src/maptel_test2.cc)

add_executable(test2 EXCLUDE_FROM_ALL ${SOURCE_TEST_2})

set_target_properties(test2 PROPERTIES OUTPUT_NAME maptel_test2)

