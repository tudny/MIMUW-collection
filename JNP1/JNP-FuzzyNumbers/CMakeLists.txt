cmake_minimum_required(VERSION 3.0)
project(Zadanie3_fuzzy)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif ()

# We don't want to see make commands.
set (CMAKE_VERBOSE_MAKEFILE ON)

# Set compiler version to GCC 10 in order for the project to work properly.
set (CMAKE_CXX_COMPILER "g++-10")

# Settings flags for compiling CPP.

set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_FLAGS "-std=c++20 -lm")

set (DEBUG_FLAGS "-g -Wall -Wextra")
set (CMAKE_CXX_FLAGS_DEBUG ${DEBUG_FLAGS})

set (RELEASE_FLAGS "-O2 -DNDEBUG")
set (CMAKE_CXX_FLAGS_RELEASE ${RELEASE_FLAGS})

# Settings source files for default executable file.

set (SOURCE_FILES
        src/fuzzy_example.cc
        src/fuzzy.cc
        src/fuzzy.h )

add_executable (${CMAKE_PROJECT_NAME} ${SOURCE_FILES})

