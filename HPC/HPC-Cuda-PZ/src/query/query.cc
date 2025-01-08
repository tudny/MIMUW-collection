#include <vector>
#include <memory>
#include <fstream>
#include <omp.h>
#include "query.h"
#include "../logger.h"
#include "worker.h"

using IndexData = std::tuple<
        std::shared_ptr<std::vector<IndexHash>>,
        std::shared_ptr<std::vector<IndexValue>>,
        char*
>;

IndexData load_index(const std::string &index_file) {
    LOGDEBUG("Loading the index_keys from the file: " << index_file);

    std::ifstream file(index_file, std::ios::binary);
    if (!file.is_open()) {
        LOGERROR("Failed to open the index_keys file: " << index_file);
        throw std::runtime_error("Failed to open the index_keys file");
    }

    uint64_t index_size;
    file.read(reinterpret_cast<char *>(&index_size), sizeof(index_size));
    LOGDEBUG("Index size: " << index_size);

    auto index_keys = std::make_shared<std::vector<IndexHash>>(index_size);
    file.read(reinterpret_cast<char *>(index_keys->data()), index_size * sizeof(IndexHash));

    auto index_values = std::make_shared<std::vector<IndexValue>>(index_size);
    file.read(reinterpret_cast<char *>(index_values->data()), index_size * sizeof(IndexValue));

    char c[2];
    file.read(c, 2);
    if (c[0] != '\n' || c[1] != '\n') {
        LOGERROR("Index file is corrupted: " << index_file);
        throw std::runtime_error("Index file is corrupted");
    }

    auto database_start_index = file.tellg();

    std::filesystem::path index_path(index_file);
    size_t index_size_bytes = std::filesystem::file_size(index_path);
    size_t data_size_bytes = index_size_bytes - database_start_index;

    char *data = new char[data_size_bytes + 1];
    file.read(data, static_cast<long>(data_size_bytes));
    for (uint64_t i = 0; i < data_size_bytes; i++) {
        if (data[i] == '\n') { data[i] = '\0'; }
    }
    data[data_size_bytes] = '\0';

    return std::make_tuple(index_keys, index_values, data);
}

bool read_query(
        const std::string &query_file,
        QueryList &queries,
        uint64_t &start,
        uint64_t max_queries = 10000000
) {
    LOGDEBUG("Reading the query file: " << query_file);

    std::ifstream file(query_file);
    if (!file.is_open()) {
        LOGERROR("Failed to open the query file: " << query_file);
        throw std::runtime_error("Failed to open the query file");
    }
    file.seekg(start);
    queries.clear();

    std::string line;
    while ((queries.size() < max_queries || max_queries == 0) && std::getline(file, line)) {
        start = file.tellg();
        if (line.empty() || line[0] == '#') {
            continue;
        }

        IndexHash query = make_clear_index();
        auto line_c = line.c_str();
        uint64_t offset = 0;
        set_index_chromosome(query, code_chromosome(line_c, offset));
        set_index_pos(query, code_position(line_c, offset));
        read_string_and_ignore(line_c, offset);
        set_index_ref(query, code_ref(line_c, offset));
        set_index_alt(query, code_alt(line_c, offset));

        queries.push_back(query);
    }

    return !queries.empty();
}

void query_index_on_single_cpu(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
) {
    for (size_t i = 0; i < queries.size(); i++) {
        const auto &query = queries[i];
        for (size_t j = 0; j < indexes->size(); j++) {
            if (is_index_same(query, indexes->at(j))) {
                results[i].offset_idx = j;
                break;
            }
        }
    }
}

void query_index_on_n_cpu(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
) {
#pragma omp parallel for
    for (size_t i = 0; i < queries.size(); i++) {
        const auto &query = queries[i];
        for (size_t j = 0; j < indexes->size(); j++) {
            if (is_index_same(query, indexes->at(j))) {
                results[i].offset_idx = j;
                break;
            }
        }
    }
}

void query_index_on_single_cpu_binary(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
) {
    for (size_t i = 0; i < queries.size(); i++) {
        const auto &query = queries[i];
        int64_t lower = -1;
        int64_t upper = indexes->size();
        while (upper - lower > 1) {
            int64_t mid = (lower + upper) / 2;
            if (query == indexes->at(mid)) {
                results[i].offset_idx = mid;
                break;
            } else if (query < indexes->at(mid)) {
                upper = mid;
            } else {
                lower = mid;
            }
        }
    }
}

void query_index_on_n_cpu_binary(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
) {
#pragma omp parallel for
    for (size_t i = 0; i < queries.size(); i++) {
        const auto &query = queries[i];
        int64_t lower = -1;
        int64_t upper = indexes->size();
        while (upper - lower > 1) {
            int64_t mid = (lower + upper) / 2;
            if (query == indexes->at(mid)) {
                results[i].offset_idx = mid;
                break;
            } else if (query < indexes->at(mid)) {
                upper = mid;
            } else {
                lower = mid;
            }
        }
    }
}

void query_index_on_single_cpu_binary_oneway(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
) {
    int64_t log = 0;
    while ((1 << log) < (int64_t) indexes->size()) {
        log++;
    }

    for (size_t i = 0; i < queries.size(); i++) {
        const auto &query = queries[i];
        int64_t begin = -1;
        for (int64_t step = 1 << log; step >= 1; step >>= 1) {
            int64_t next = begin + step;
            if ((next < (int) indexes->size()) && (query >= indexes->at(next))) {
                begin = next;
            }
        }
        if (begin != -1) {
            if (query == indexes->at(begin)) {
                results[i].offset_idx = begin;
            }
        }
    }
}

void query_index_on_n_cpu_binary_oneway(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
) {
    int64_t log = 0;
    while ((1 << log) < (int64_t) indexes->size()) {
        log++;
    }

#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t) queries.size(); i++) {
        const auto &query = queries[i];
        int64_t begin = -1;
        for (int64_t step = 1 << log; step >= 1; step >>= 1) {
            int64_t next = begin + step;
            if ((next < (int64_t) indexes->size()) && (query >= indexes->at(next))) {
                begin = next;
            }
        }
        if (begin != -1) {
            if (query == indexes->at(begin)) {
                results[i].offset_idx = begin;
            }
        }
    }
}

void unknown_query_resolution(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
) {
    (void) indexes;
    (void) queries;
    (void) results;
    throw std::runtime_error("Unknown query resolution strategy");
}

void resolve_query(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results,
        const QueryStrategy &strategy
) {
    query_resolution resolution = unknown_query_resolution;
    switch (strategy) {
        case CPU_1CORE_LINEAR:
            resolution = query_index_on_single_cpu;
            break;
        case CPU_1CORE_BINARY:
            resolution = query_index_on_single_cpu_binary;
            break;
        case CPU_1CORE_BINARY_ONEWAY:
            resolution = query_index_on_single_cpu_binary_oneway;
            break;
        case CPU_NCORE_LINEAR:
            resolution = query_index_on_n_cpu;
            break;
        case CPU_NCORE_BINARY:
            resolution = query_index_on_n_cpu_binary;
            break;
        case CPU_NCORE_BINARY_ONEWAY:
            resolution = query_index_on_n_cpu_binary_oneway;
            break;
        case GPU_LINEAR:
            resolution = query_index_on_gpu_linear;
            break;
        case GPU_BINARY:
            resolution = query_index_on_gpu_binary;
            break;
        case GPU_BINARY_ONEWAY:
            resolution = query_index_on_gpu_binary_oneway;
            break;
        case BENCHMARK_ALL:
            goto benchmark_all;
    }

    TIMEIT(query_strategy_name(strategy), resolution(indexes, queries, results));
    return;

    benchmark_all:
    for (const auto &[name, queryStrategy]: STRATEGY_MAP) {
        if (queryStrategy == BENCHMARK_ALL) {
            continue;
        }
        results.clear();
        resolve_query(indexes, queries, results, queryStrategy);
    }
}

void query(const QueryParams &params) {
    LOGDEBUG("Querying the index: " << "Query: " << params.query
                                    << " Index: " << params.index
                                    << " Output: " << params.output
                                    << " Strategy: " << query_strategy_name(params.strategy));

    validate_file_exists(params.index);
    validate_file_exists(params.query);
    DO_IN_DEBUG(dump_file_info(params.index));
    DO_IN_DEBUG(dump_file_info(params.query));

    auto [index_keys, index_values, buffer] = load_index(params.index);

    LOGDEBUG("Loaded the index: " << "Index size: " << index_keys->size());

    std::vector<IndexHash> queries;
    uint64_t start = 0;

    create_path_if_not_exists(params.output);
    std::ofstream out_file(params.output, std::ios::binary);

    while (read_query(params.query, queries, start, 0)) {
        LOGDEBUG("Read the queries: " << "Number of queries: " << queries.size());

        std::vector<QueryResult> results{queries.size()};
        resolve_query(index_keys, queries, results, params.strategy);

        LOGDEBUG("Queried the index: " << "Number of results: " << results.size());

        DO_IN_DEBUG(
                int64_t non_empty_results = 0;
                for (const auto &result: results) {
                    if (result.offset_idx != NO_RESULT) {
                        non_empty_results++;
                    }
                }
                LOGDEBUG("Number of non-empty results: " << non_empty_results);
        )

        for (const auto &result: results) {
            if (result.offset_idx == NO_RESULT) {
                continue;
            }
            auto offset = index_values->at(result.offset_idx);
            auto line = buffer + offset;
            out_file << line << '\n';
        }
    }

    free(buffer);
}
