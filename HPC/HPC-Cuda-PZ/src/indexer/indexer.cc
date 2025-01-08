#include <fstream>
#include <vector>
#include <cassert>
#include <queue>
#include "indexer.h"
#include "../logger.h"

uint64_t count_entries(const std::string &file) {
    std::ifstream in_file(file);
    uint64_t count = 0;
    std::string line;
    while (std::getline(in_file, line)) {
        if (!line.empty() && line[0] != '#') {
            ++count;
        }
    }
    in_file.close();
    return count;
}

using Index = std::pair<IndexHash, IndexValue>;

std::string index_hash_to_string(const IndexHash &index_hash, const IndexValue &index_value) {
    return "chr=" + std::to_string(get_index_chromosome(index_hash)) + " "
           + "pos=" + std::to_string(get_index_pos(index_hash)) + " "
           + "ref=" + std::to_string(get_index_ref(index_hash)) + " "
           + "alt=" + std::to_string(get_index_alt(index_hash)) + " "
           + "offset=" + std::to_string(index_value);
}

using UnsortedLocalIndexQueue = std::priority_queue<Index, std::vector<Index>, std::greater<>>;

void add_remaining_unsorted_to_sorted(
        UnsortedLocalIndexQueue &unsorted_indexes,
        std::vector<IndexHash> &sorted_indexes_keys,
        std::vector<IndexValue> &sorted_indexes_values
) {
    while (!unsorted_indexes.empty()) {
        auto next_index = unsorted_indexes.top();
        sorted_indexes_keys.push_back(next_index.first);
        sorted_indexes_values.push_back(next_index.second);
        unsorted_indexes.pop();
    }
}

void add_index_to_unsorted_and_update_sorted(
        UnsortedLocalIndexQueue &unsorted_indexes,
        std::vector<IndexHash> &sorted_indexes_keys,
        std::vector<IndexValue> &sorted_indexes_values,
        const Index &index
) {
    if (unsorted_indexes.empty() || !is_index_same_semigroup(unsorted_indexes.top().first, index.first)) {
        add_remaining_unsorted_to_sorted(unsorted_indexes, sorted_indexes_keys, sorted_indexes_values);
    }
    unsorted_indexes.push(index);
}

void index_file(const IndexingParams &params) {
    uint64_t estimated_entries = count_entries(params.database);
    LOGDEBUG("Estimated entries: " << estimated_entries);

    std::string line;
    std::ifstream in_file(params.database);
    std::vector<IndexHash> indexes_keys;
    std::vector<IndexValue> indexes_values;

    // This shouldn't be exceeded 2GB even on the whole database.
    // We can safely reserve this amount of memory.
    indexes_keys.reserve(estimated_entries);
    indexes_values.reserve(estimated_entries);

    // We use a priority queue to keep the indexes sorted in the same semigroup.
    // There are at most 16 elements with the same semigroup.
    // This is a small number, and we can afford to keep them in memory.
    std::vector<Index> _temp;
    _temp.reserve(16);
    UnsortedLocalIndexQueue unsorted_indexes{std::greater<>(), std::move(_temp)};

    while (true) {
        uint64_t offset = 0;
        Index index = {make_clear_index(), make_clear_index_value()};
        index.second = in_file.tellg();
        if (!std::getline(in_file, line)) {
            break;
        }

        if (line.empty() || line[0] == '#') {
            continue;
        }

        auto line_c = line.c_str();
        set_index_chromosome(index.first, code_chromosome(line_c, offset));
        set_index_pos(index.first, code_position(line_c, offset));
        set_index_ref(index.first, code_ref(line_c, offset));
        set_index_alt(index.first, code_alt(line_c, offset));

        add_index_to_unsorted_and_update_sorted(unsorted_indexes, indexes_keys, indexes_values, index);
    }
    add_remaining_unsorted_to_sorted(unsorted_indexes, indexes_keys, indexes_values);

    create_path_if_not_exists(params.index);
    std::ofstream index_file(params.index, std::ios::binary | std::ios::out);
    uint64_t number_of_indexes = indexes_keys.size();
    index_file.write(reinterpret_cast<const char *>(&number_of_indexes), sizeof(number_of_indexes));
    index_file.write(reinterpret_cast<const char *>(indexes_keys.data()), number_of_indexes * sizeof(IndexHash));
    index_file.write(reinterpret_cast<const char *>(indexes_values.data()), number_of_indexes * sizeof(IndexValue));

    LOGDEBUG("Indexing done: " << "Number of indexes: " << number_of_indexes);

    assert(number_of_indexes == estimated_entries);
    DO_IN_DEBUG(
            for (uint64_t idx = 0; idx + 1 < number_of_indexes; ++idx) {
                if (!(indexes_keys[idx] <= indexes_keys[idx + 1])) {
                    LOGERROR(
                            "Index is not sorted: idx=" << idx
                                                        << " " << index_hash_to_string(indexes_keys[idx],
                                                                                       indexes_values[idx])
                                                        << " " << index_hash_to_string(indexes_keys[idx + 1],
                                                                                       indexes_values[idx + 1])
                    );
                }
                assert(indexes_keys[idx] <= indexes_keys[idx + 1]);
            }
    );

    in_file.close();
    // For sanity, we add two newlines to separate the indexes from the database.
    index_file.write("\n\n", 2);

    in_file.open(params.database);
    index_file << in_file.rdbuf();
    in_file.close();
    index_file.close();
}

void index(const IndexingParams &params) {
    LOGDEBUG("Indexing the database: " << "Database: " << params.database
                                       << " Index: " << params.index);

    validate_file_exists(params.database);
    DO_IN_DEBUG(dump_file_info(params.database));
    index_file(params);
}
