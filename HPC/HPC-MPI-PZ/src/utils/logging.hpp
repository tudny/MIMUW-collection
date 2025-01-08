#pragma once

#include <chrono>
#include <cstring>
#include <ctime>
#include <iostream>
#include <mpi.h>

// Logging is disabled in release mode as it can have a significant impact on performance
#ifndef NDEBUG
#define LOG(msg, status)                                                       \
  do {                                                                         \
    int my_super_mega_giga_rank;                                                      \
    MPI_Comm_rank(MPI_COMM_WORLD, &my_super_mega_giga_rank);                          \
    if (my_super_mega_giga_rank == 0) {                                                                           \
    auto current_time = std::chrono::system_clock::now();                      \
    std::string path{__FILE__};                                                \
    path = path.substr(path.find_last_of("/\\") + 1);                          \
    std::time_t current_time_t =                                               \
        std::chrono::system_clock::to_time_t(current_time);                    \
    std::cerr << std::strtok(std::ctime(&current_time_t), "\n")                \
              << " [" << status << "] "                                        \
              << path << ":" << __LINE__ << " - " << __func__ << " - "     \
              << msg << std::endl;       }                                      \
  } while (false)
#else
#define LOG(msg, status)
#endif

#define LOGINFO(msg) LOG(msg, "INFO ")
#define LOGERROR(msg) LOG(msg, "ERROR")
#define LOGWARN(msg) LOG(msg, "WARN ")
#define LOGDEBUG(msg) LOG(msg, "DEBUG")

#define LOGTIME(tag) \
    do { \
        int my_super_mega_giga_rank; \
        MPI_Comm_rank(MPI_COMM_WORLD, &my_super_mega_giga_rank); \
        if (my_super_mega_giga_rank == 0) {                      \
            auto current_time = std::chrono::system_clock::now(); \
            auto duration = current_time.time_since_epoch(); \
            auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count(); \
            std::cerr << millis << " [TIME] " << tag << std::endl; \
        } \
    } while (false)
