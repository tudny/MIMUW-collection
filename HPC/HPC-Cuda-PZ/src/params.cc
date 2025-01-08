#include "params.h"
#include "logger.h"

#define CHECK_ARGC(condition)                       \
do {                                                \
    if (!(condition)) {                             \
        LOGDEBUG("Invalid number of arguments");    \
        return std::monostate();                    \
    }                                               \
} while (false)

InputParams parse_args(int argc, char **argv) {
    LOGDEBUG("Parsing the input arguments");

    CHECK_ARGC(argc == 4 || argc == 5);

    if (std::string(argv[1]) == "-i") {
        CHECK_ARGC(argc == 4);
        LOGDEBUG("Parsing the input arguments as IndexingParams");
        IndexingParams params;
        params.database = argv[2];
        params.index = argv[3];
        return params;
    } else {
        CHECK_ARGC(argc == 5 || argc == 4);
        LOGDEBUG("Parsing the input arguments as QueryParams");
        QueryParams params;
        params.query = argv[1];
        params.index = argv[2];
        params.output = argv[3];
        if (argc == 5) {
            params.strategy = parse_query_strategy(argv[4]);
        }
        return params;
    }
}

QueryStrategy parse_query_strategy(const std::string &strategy_name) {
    LOGDEBUG("Parsing the query strategy: " << strategy_name);

    auto it = STRATEGY_MAP.find(strategy_name);
    if (it == STRATEGY_MAP.end()) {
        LOGDEBUG("Unknown query strategy: " << strategy_name);
        return DEFAULT_QUERY_STRATEGY;
    }

    return it->second;
}

std::string query_strategy_name(QueryStrategy strategy) {
    auto it = STRATEGY_NAMES.find(strategy);
    if (it == STRATEGY_NAMES.end()) {
        return "UNKNOWN";
    }

    return it->second;
}
