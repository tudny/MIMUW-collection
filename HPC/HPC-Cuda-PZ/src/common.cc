#include <filesystem>
#include <cmath>
#include <cassert>
#include "common.h"
#include "error.h"

// source: https://en.cppreference.com/w/cpp/filesystem/file_size
std::ostream &operator<<(std::ostream &os, HumanReadable hr) {
    int o{};
    auto mantissa = (double) hr.size;
    for (; mantissa >= 1024.; mantissa /= 1024., ++o);
    os << std::ceil(mantissa * 10.) / 10. << "BKMGTPE"[o];
    return o ? os << "B (" << hr.size << ')' : os;
}

void validate_file_exists(const std::string &file) {
    LOGDEBUG("Validating file exists: " << "File: " << file);

    if (!fs::exists(file)) {
        PRINT_ERROR("File does not exist: " + file);
    }

    LOGDEBUG("File exists: " << file);
}

void dump_file_info(const std::string &file) {
    LOGDEBUG("-- Dumping file information: " << "File: " << file);

    fs::path path(file);
    LOGDEBUG("File name: " << path.filename());
    LOGDEBUG("File extension: " << path.extension());
    LOGDEBUG("File parent path: " << path.parent_path());
    LOGDEBUG("File root path: " << path.root_path());
    LOGDEBUG("File relative path: " << path.relative_path());
    LOGDEBUG("File stem: " << path.stem());
    LOGDEBUG("File size: " << HumanReadable{fs::file_size(file)});
}

bool is_end_of_data(const char *buff, uint64_t offset) {
    return buff[offset] == '\t' || buff[offset] == '\n' || buff[offset] == '\r' || buff[offset] == '\0';
}

void move_offset(const char *buff, uint64_t &offset) {
    (void) buff;
    DO_IN_DEBUG(assert(is_end_of_data(buff, offset)));
    ++offset;
}

uint8_t code_chromosome(const char *buff, uint64_t &offset) {
    uint8_t result = 0;
    if (buff[offset] == 'X') {
        result = 23;
        ++offset;
    } else if (buff[offset] == 'Y') {
        result = 24;
        ++offset;
    } else if (buff[offset] == 'M') {
        result = 25;
        ++offset;
    } else {
        result = code_uint<uint8_t>(buff, offset);
    }
    return result;
}

uint64_t code_position(const char *buff, uint64_t &offset) {
    return code_uint<uint64_t>(buff, offset);
}

bool is_actg(char c) {
    return c == 'A' || c == 'C' || c == 'T' || c == 'G';
}

uint8_t code_actg(const char *buff, uint64_t &offset) {
    uint8_t result_ascii = buff[offset++];
    DO_IN_DEBUG(assert(is_actg(result_ascii)));
    move_offset(buff, offset);
    switch (result_ascii) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'T':
            return 2;
        case 'G':
            return 3;
        default:
            PRINT_ERROR("Invalid nucleotide: " + std::to_string(result_ascii));
    }
}

uint8_t code_ref(const char *buff, uint64_t &offset) {
    return code_actg(buff, offset);
}

uint8_t code_alt(const char *buff, uint64_t &offset) {
    return code_actg(buff, offset);
}

void read_string_and_ignore(const char *buff, uint64_t &offset) {
    while (!is_end_of_data(buff, offset)) {
        ++offset;
    }
    move_offset(buff, offset);
}

IndexHash make_clear_index() {
    return 0;
}

IndexHash make_index(uint8_t chromosome, uint64_t pos, uint8_t ref, uint8_t alt) {
    IndexHash index = make_clear_index();
    set_index_chromosome(index, chromosome);
    set_index_pos(index, pos);
    set_index_ref(index, ref);
    set_index_alt(index, alt);
    return index;
}

struct mask_t {
    uint64_t bits;
    uint64_t offset;
    uint64_t mask;

    void constexpr off(uint64_t off) {
        offset += off;
        mask <<= off;
    }
};

static constexpr mask_t make_mask_pair(uint64_t bits) {
    return {bits, 0, ((1ULL << bits) - 1ULL)};
}

static constexpr mask_t make_mask_pair(uint64_t bits, mask_t prev_mask) {
    uint64_t offset = prev_mask.bits + prev_mask.offset;
    auto mask = make_mask_pair(bits);
    mask.off(offset);
    return mask;
}

static constexpr mask_t ALT_MASK = make_mask_pair(2);
static constexpr mask_t REF_MASK = make_mask_pair(2, ALT_MASK);
static constexpr mask_t POS_MASK = make_mask_pair(55, REF_MASK);
static constexpr mask_t CHR_MASK = make_mask_pair(5, POS_MASK);

static_assert((POS_MASK.bits + CHR_MASK.bits + REF_MASK.bits + ALT_MASK.bits) == 64);
static_assert((POS_MASK.mask & CHR_MASK.mask & REF_MASK.mask & ALT_MASK.mask) == 0);
static_assert((POS_MASK.mask | CHR_MASK.mask | REF_MASK.mask | ALT_MASK.mask) == 0xFFFFFFFFFFFFFFFF);

static void _modify_index(IndexHash &index, mask_t mask, uint64_t value) {
    index = (index & ~mask.mask) | ((uint64_t) value << mask.offset);
}

static uint64_t _get_index(IndexHash index, mask_t mask) {
    return (index & mask.mask) >> mask.offset;
}

void set_index_chromosome(IndexHash &index, uint8_t chromosome) {
    _modify_index(index, CHR_MASK, chromosome);
}

void set_index_pos(IndexHash &index, uint64_t pos) {
    _modify_index(index, POS_MASK, pos);
}

void set_index_ref(IndexHash &index, uint8_t ref) {
    _modify_index(index, REF_MASK, ref);
}

void set_index_alt(IndexHash &index, uint8_t alt) {
    _modify_index(index, ALT_MASK, alt);
}

uint8_t get_index_chromosome(IndexHash index) {
    return _get_index(index, CHR_MASK);
}

uint64_t get_index_pos(IndexHash index) {
    return _get_index(index, POS_MASK);
}

uint8_t get_index_ref(IndexHash index) {
    return _get_index(index, REF_MASK);
}

uint8_t get_index_alt(IndexHash index) {
    return _get_index(index, ALT_MASK);
}

IndexValue make_clear_index_value() {
    return 0;
}

IndexValue make_index_offset(uint64_t offset) {
    return offset;
}

void set_index_offset(IndexValue &index, uint64_t offset) {
    index = offset;
}

uint64_t get_index_offset(IndexValue index) {
    return index;
}

// IndexHash have the same semigroup if they have the same chromosome and position
bool is_index_same_semigroup(IndexHash a, IndexHash b) {
    uint64_t mask = CHR_MASK.mask | POS_MASK.mask;
    return (a & mask) == (b & mask);
}

void create_path_if_not_exists(const std::string &file_path) {
    LOGDEBUG("Creating path if not exists: " << "Path: " << file_path);

    fs::path path(file_path);
    if (path.has_parent_path()) {
        fs::create_directories(path.parent_path());
    }

    LOGDEBUG("Path created: " << file_path);
}

template<typename T>
T code_uint(const char *buff, uint64_t &offset) {
    T result = 0;
    while (!is_end_of_data(buff, offset)) {
        result = result * 10 + (buff[offset] - '0');
        ++offset;
    }
    move_offset(buff, offset);
    return result;
}
