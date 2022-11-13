cmake_minimum_required(VERSION 3.0)
project(Zadanie6_rover)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif ()

# We don't want to see make commands.
set (CMAKE_VERBOSE_MAKEFILE ON)

# Set compiler version to GCC 11 in order for the project to work properly.
set (CMAKE_CXX_COMPILER "g++-11")

# Settings flags for compiling CPP.

set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_FLAGS "-std=c++20")

set (DEBUG_FLAGS "-g -Wall -Wextra")
set (CMAKE_CXX_FLAGS_DEBUG ${DEBUG_FLAGS})

set (RELEASE_FLAGS "-O2")
set (CMAKE_CXX_FLAGS_RELEASE ${RELEASE_FLAGS})

# Adding the rover library
add_library(Rover
        src/rover.cc
        src/command.cc)

# Settings source files for default executable file.
add_executable(${CMAKE_PROJECT_NAME} src/rover_example.cc)
target_link_libraries(${CMAKE_PROJECT_NAME} Rover)

# Setting source files for tests.
set(TEST_NAME ${CMAKE_PROJECT_NAME}_test)
add_executable(${TEST_NAME} src/rover_test.cc)
target_link_libraries(${TEST_NAME} Rover)
set_target_properties(${TEST_NAME} PROPERTIES OUTPUT_NAME ${TEST_NAME})
