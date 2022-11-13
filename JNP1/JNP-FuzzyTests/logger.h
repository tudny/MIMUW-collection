#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <vector>

namespace Logger {

    inline void log(const std::string &message) {
        std::cout << "> " << message << std::endl;
    }

    template<typename... Messages>
    inline void log(const std::string &message, Messages... messages) {
        log(message);
        log(messages...);
    }

    inline void log(const std::vector<std::string> &messages) {
        for (const auto &message : messages)
            log(message);
    }

}

#endif // LOGGER_H
