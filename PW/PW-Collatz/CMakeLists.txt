cmake_minimum_required(VERSION 3.0)
project(PW-Collatz)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif ()

# We want to see make commands.
set (CMAKE_VERBOSE_MAKEFILE ON)

set (CMAKE_CXX_STANDARD 14)

add_executable(main main.cpp teams.cpp)
add_executable(new_process new_process.cpp)

target_link_libraries(main pthread)
target_link_libraries(new_process pthread)

install(TARGETS DESTINATION .)
