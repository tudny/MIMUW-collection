cmake_minimum_required(VERSION 3.10)
project(HPC_MPI_PZ)

if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
    cmake_policy(SET CMP0135 NEW)
endif ()

set(EXEC_NAME matmul)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif ()

cmake_host_system_information(RESULT CMAKE_HOSTNAME QUERY HOSTNAME)

if (CMAKE_HOSTNAME MATCHES "okeanos.*")
    message(STATUS "Running on okeanos, setting Cray compiler")
    set(CMAKE_C_COMPILER /opt/cray/pe/craype/2.6.1/bin/cc)
    set(CMAKE_CXX_COMPILER /opt/cray/pe/craype/2.6.1/bin/CC)
else ()
    message(STATUS "Running on local machine, leaving default compiler")
endif ()

set(CMAKE_VERBOSE_MAKEFILE OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wconversion -std=c++17 -Wno-cast-function-type")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

set(SOURCE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/src)
file(GLOB_RECURSE SOURCE_FILES CONFIGURE_DEPENDS ${SOURCE_PATH}/*.cpp ${SOURCE_PATH}/*.hpp)
add_executable(${EXEC_NAME} ${SOURCE_FILES})

find_package(MPI QUIET)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

target_link_libraries(${EXEC_NAME} LINK_PUBLIC ${MPI_LIBRARIES} stdc++fs)
include_directories(${MPI_INCLUDE_PATH})

if (CMAKE_HOSTNAME MATCHES "okeanos.*")
    message("Omitting test compilation on Okeanos")
    return()
endif ()

# For Unit Tests
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

set(TEST_PATH ${CMAKE_CURRENT_SOURCE_DIR}/test)
set(TEST_EXEC_NAME ${PROJECT_NAME}_test)
file(GLOB_RECURSE TEST_FILES CONFIGURE_DEPENDS ${TEST_PATH}/*.cpp ${TEST_PATH}/*.hpp*)
set(SOURCE_FILES_WITHOUT_MAIN ${SOURCE_FILES})
list(FILTER SOURCE_FILES_WITHOUT_MAIN EXCLUDE REGEX ".*matmul.cpp$")
add_executable(
        ${TEST_EXEC_NAME}
        ${TEST_FILES}
        ${SOURCE_FILES_WITHOUT_MAIN}
)
target_link_libraries(
        ${TEST_EXEC_NAME}
        GTest::gtest_main
        GTest::gmock_main
        ${MPI_LIBRARIES}
)

include(GoogleTest)
gtest_discover_tests(${TEST_EXEC_NAME})

set_target_properties(${TEST_EXEC_NAME} PROPERTIES COMPILE_FLAGS "-DUNIT_TEST_DISABLE_MPI")

configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/test/assets.hpp.in
        ${CMAKE_CURRENT_BINARY_DIR}/test/assets.hpp
)
include_directories(${TEST_EXEC_NAME}
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test
        PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/test
)
