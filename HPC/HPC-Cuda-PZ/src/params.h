#ifndef GPUGENV_PARAMS_H
#define GPUGENV_PARAMS_H

#include <string>
#include <variant>
#include <map>

template<class... Ts>
struct overloaded : Ts ... {
    using Ts::operator()...;
};
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

struct IndexingParams {
    std::string database;
    std::string index;
};

#define STRATEGIES \
    X(CPU_1CORE_LINEAR) \
    X(CPU_NCORE_LINEAR) \
    X(CPU_1CORE_BINARY) \
    X(CPU_NCORE_BINARY) \
    X(CPU_1CORE_BINARY_ONEWAY) \
    X(CPU_NCORE_BINARY_ONEWAY) \
    X(GPU_LINEAR) \
    X(GPU_BINARY)  \
    X(GPU_BINARY_ONEWAY)  \
    X(BENCHMARK_ALL)

#define X(name) name,
enum QueryStrategy {
    STRATEGIES
};
#undef X

#define X(name) {#name, name},
static const std::map<std::string, QueryStrategy> STRATEGY_MAP = {
        STRATEGIES
};
#undef X

#define X(name) {name, #name},
static const std::map<QueryStrategy, std::string> STRATEGY_NAMES = {
        STRATEGIES
};
#undef X

static const QueryStrategy DEFAULT_QUERY_STRATEGY = GPU_BINARY_ONEWAY;

QueryStrategy parse_query_strategy(const std::string &strategy_name);

std::string query_strategy_name(QueryStrategy strategy);

struct QueryParams {
    std::string query;
    std::string index;
    std::string output;
    QueryStrategy strategy = DEFAULT_QUERY_STRATEGY;
};

using InputParams = std::variant<IndexingParams, QueryParams, std::monostate>;

InputParams parse_args(int argc, char **argv);

#endif //GPUGENV_PARAMS_H
