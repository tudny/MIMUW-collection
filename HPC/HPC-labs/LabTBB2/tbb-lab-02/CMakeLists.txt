cmake_minimum_required (VERSION 3.0)
project (MIMUW-tbb-lab1 CXX)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")

find_package ( TBB REQUIRED )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-g -Wall")

add_executable(pardo pardo.cpp)
target_link_libraries(pardo tbb_debug)

add_executable(pardo_th_sp pardo_th_sp.cpp)
target_link_libraries(pardo_th_sp tbb_debug)

add_executable(dij dij.cpp)
target_link_libraries(dij tbb_debug)
