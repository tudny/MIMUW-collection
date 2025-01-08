#ifndef GPUGENV_LOGGER_H
#define GPUGENV_LOGGER_H

#include <chrono>
#include <cstring>
#include <ctime>
#include <iostream>

#define LOG(msg, status)                                                     \
  do {                                                                         \
    auto current_time = std::chrono::system_clock::now();                      \
    std::time_t current_time_t =                                               \
        std::chrono::system_clock::to_time_t(current_time);                    \
    std::cerr << std::strtok(std::ctime(&current_time_t), "\n")                \
              << " [" << status << "] "                                        \
              << __FILE__ << ":" << __LINE__ << " - " << __func__ << " - "     \
              << msg << std::endl;                                             \
  } while (false)

#ifdef NDEBUG
#define LOGDEBUG(msg)                                                          \
  do {                                                                         \
  } while (false)
#else
#define LOGDEBUG(msg) LOG(msg, "DEBUG")
#endif

#define LOGINFO(msg) LOG(msg, "INFO")
#define LOGERROR(msg) LOG(msg, "ERROR")

#endif // GPUGENV_LOGGER_H
