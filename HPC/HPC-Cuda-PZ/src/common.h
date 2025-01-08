#ifndef GPUGENV_COMMON_H
#define GPUGENV_COMMON_H

#include <filesystem>

#ifndef NDEBUG
#define DO_IN_DEBUG(code) code
#else
#define DO_IN_DEBUG(code)
#endif

namespace fs = std::filesystem;

struct HumanReadable {
    std::uintmax_t size{};

private:
    friend std::ostream &operator<<(std::ostream &os, HumanReadable hr);
};

void validate_file_exists(const std::string &file);

void create_path_if_not_exists(const std::string &file_path);

void dump_file_info(const std::string &file);

bool is_end_of_data(const char *buff, uint64_t offset);

using IndexHash = uint64_t;

IndexHash make_clear_index();

IndexHash make_index(uint8_t chromosome, uint64_t pos, uint8_t ref, uint8_t alt);

void set_index_chromosome(IndexHash &index, uint8_t chromosome);

void set_index_pos(IndexHash &index, uint64_t pos);

void set_index_ref(IndexHash &index, uint8_t ref);

void set_index_alt(IndexHash &index, uint8_t alt);

uint8_t get_index_chromosome(IndexHash index);

uint64_t get_index_pos(IndexHash index);

uint8_t get_index_ref(IndexHash index);

uint8_t get_index_alt(IndexHash index);

bool is_index_same_semigroup(IndexHash a, IndexHash b);

inline bool is_index_same(IndexHash a, IndexHash b) {
    return a == b;
}

// IndexValue encodes data in the following way:
// first 64 bits are for offset
using IndexValue = uint64_t;

IndexValue make_clear_index_value();

IndexValue make_index_offset(uint64_t offset);

void set_index_offset(IndexValue &index, uint64_t offset);

uint64_t get_index_offset(IndexValue index);

template<typename T>
T code_uint(const char *buff, uint64_t &offset);

void read_string_and_ignore(const char *buff, uint64_t &offset);

// 1-22 -> 1-22, X -> 23, Y -> 24, M -> 25
uint8_t code_chromosome(const char *buff, uint64_t &offset);

uint64_t code_position(const char *buff, uint64_t &offset);

bool is_actg(char c);

uint8_t code_actg(const char *buff, uint64_t &offset);

uint8_t code_ref(const char *buff, uint64_t &offset);

uint8_t code_alt(const char *buff, uint64_t &offset);

#define TIMEIT_BEGIN(start) auto start = std::chrono::high_resolution_clock::now()
#define TIMEIT_END(start) do { \
    auto end = std::chrono::high_resolution_clock::now(); \
    LOGINFO("[TIMEIT] " << #start << " - Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms"); \
} while (false)

#define TIMEIT(mss, code) do { \
    auto start = std::chrono::high_resolution_clock::now(); \
    code; \
    auto end = std::chrono::high_resolution_clock::now(); \
    LOGINFO("[TIMEIT] " << mss << " - Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms"); \
} while (false)

#endif // GPUGENV_COMMON_H
