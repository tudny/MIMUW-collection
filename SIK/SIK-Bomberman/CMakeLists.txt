cmake_minimum_required(VERSION 3.16)
project(SIK-Bomberman)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif ()

# We want to see make commands.
set (CMAKE_VERBOSE_MAKEFILE ON)

cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
message("Hostname is: " ${_host_name})

# Settings GCC 11.2 compiler
if (${_host_name} STREQUAL "students")
    message ("Setting STUDENTS gcc11.2 path")
    message ("Before running your program run: export LD_LIBRARY_PATH=/opt/gcc-11.2/lib64")
    set (LD_LIBRARY_PATH "/opt/gcc-11.2/lib64")
    set (CMAKE_CXX_COMPILER "/opt/gcc-11.2/bin/g++-11.2")
else ()
    message (WARNING
            "You are trying to use local g++11.2 compiler.\
            Before doing so, change compiler PATH and LD_LIBRARY_PATH accordingly\
            or CMAKE this project on STUDENTS machine.")

    message ("Setting LOCAL gcc11.2 path")
    message ("Before running your program run: export LD_LIBRARY_PATH=/home/tudny/.GCC-11.2/lib64")
    set (LD_LIBRARY_PATH "/home/tudny/.GCC-11.2/lib64")
    set (CMAKE_CXX_COMPILER "/home/tudny/.GCC-11.2/bin/g++")
endif()

# Set flags
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_FLAGS "-Wall -Wextra -Wconversion -Werror -std=gnu++20 -pthread -DBOOST_LOG_DYN_LINK")
set (CMAKE_CXX_FLAGS_DEBUG "-g")
set (CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

# Add executable client
set (CLIENT_EXECUTABLE robots-client)
add_executable (${CLIENT_EXECUTABLE}
        src-client/main.cpp
        src-client/utils/arguments.cpp
        src-client/utils/arguments.h
        src-client/serializer/serializer.cpp
        src-client/serializer/serializer.h
        src-client/gui/Connection.cpp
        src-client/gui/Connection.h
        src-client/server/Connection.cpp
        src-client/server/Connection.h
        src-client/bridge/bridge-gui-server/GuiServerBridge.cpp
        src-client/bridge/bridge-gui-server/GuiServerBridge.h
        src-client/utils/utils.h
        src-client/bridge/bridge-server-gui/ServerGuiBridge.cpp
        src-client/bridge/bridge-server-gui/ServerGuiBridge.h
        src-client/bridge/Bridge.cpp
        src-client/bridge/Bridge.h
        src-client/game/GameState.cpp
        src-client/game/GameState.h)

## Add executable server
set (SERVER_EXECUTABLE robots-server)
add_executable (${SERVER_EXECUTABLE}
        src-server/main.cpp
        src-server/utils/arguments.cpp
        src-server/utils/arguments.h
        src-server/utils/RNG.cpp
        src-server/utils/RNG.h
        src-server/client/Connection.cpp
        src-server/client/Connection.h
        src-server/client/Acceptor.cpp
        src-server/client/Acceptor.h
        src-server/threads/set.h
        src-server/client/Sender.cpp
        src-server/client/Sender.h
        src-server/client/Receiver.cpp
        src-server/client/Receiver.h
        src-server/game/GameLoop.cpp
        src-server/game/GameLoop.h
        src-server/game/models/Player.h
        src-server/game/models/Position.h
        src-server/game/models/Networkable.h
        src-server/game/models/Numbers.h
        src-server/game/models/String.h
        src-server/game/models/models.h
        src-server/game/models/Bomb.h
        src-server/game/models/Containers.h
        src-server/game/models/ClientMessage.h
        src-server/game/GameState.cpp
        src-server/game/GameState.h
        src-server/game/models/Direction.h
        src-server/threads/queue.h
        src-server/client/PlayerInstance.h
        src-server/game/models/ServerMessage.h)

# Setup BOOST library
find_package(Boost 1.74 COMPONENTS program_options system thread log log_setup REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost Not found")
endif()

target_link_libraries(${CLIENT_EXECUTABLE} LINK_PUBLIC ${Boost_LIBRARIES})
target_link_libraries(${SERVER_EXECUTABLE} LINK_PUBLIC ${Boost_LIBRARIES})
