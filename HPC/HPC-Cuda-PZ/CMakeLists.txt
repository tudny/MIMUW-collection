cmake_minimum_required(VERSION 3.22)
if (NOT DEFINED CMAKE_CUDA_COMPILER)
    set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
endif()
project(gpugenv LANGUAGES CXX CUDA)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fopenmp -std=c++17")
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3 -DNDEBUG -Xptxas -O3")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -Wall,-Wextra,-fopenmp")
set(CMAKE_CUDA_FLAGS_DEBUG_INIT "-g")
set(CMAKE_CUDA_FLAGS_RELEASE_INIT "-O3 -DNDEBUG -Xptxas -O3")

add_executable(
        ${CMAKE_PROJECT_NAME}
        src/main.cc
        src/logger.h
        src/params.cc
        src/params.h
        src/indexer/indexer.cc
        src/indexer/indexer.h
        src/query/query.cc
        src/query/query.h
        src/error.h
        src/error.cc
        src/common.h
        src/common.cc
        src/query/worker.h
        src/query/worker.cu
)
