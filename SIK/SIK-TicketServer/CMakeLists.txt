cmake_minimum_required(VERSION 3.21)
project(UDPTicketServer)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif ()

# We do want to see make commands.
set (CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_FLAGS "-Wall -Wextra -Wno-implicit-fallthrough")

set(CMAKE_C_FLAGS_DEBUG "-g -D SAFE_EXIT")
set(CMAKE_C_FLAGS_RELEASE "-O2")

set(SOURCE_FILES ticket_server.c)

add_executable(ticket_server ${SOURCE_FILES})
