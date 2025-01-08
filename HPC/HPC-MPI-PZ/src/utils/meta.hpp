#pragma once

#include <cstdint>
#include <optional>
#include <tuple>

using ProcessIdentifier = std::tuple<uint32_t, uint32_t, uint32_t>;
using IndexRange = std::pair<uint32_t, uint32_t>;

extern ProcessIdentifier ProcessIdentifierMASTER;

enum SplitOverMatrix {
    A, B, I_DONT_CARE
};

/// p_r * p_c * l = number_of_processes
/// also p_r = p_c
/// if layers is not provided we can assume l = 1
struct Topology {
    /// Processes rows
    uint32_t p_r;
    /// Processes columns
    uint32_t p_c;
    /// Processes layers
    uint32_t l;

    Topology(uint32_t number_of_processes, std::optional<uint32_t> layers);

    void validate_matrix_size(uint32_t rows, uint32_t cols) const;

    [[nodiscard]] ProcessIdentifier get_process_identifier(uint32_t process_rank) const;

    [[nodiscard]] ProcessIdentifier get_process_identifier_int(int process_rank) const;

    [[nodiscard]] uint32_t get_process_rank(const ProcessIdentifier &process_identifier) const;

    [[nodiscard]] int get_process_rank_int(const ProcessIdentifier &process_identifier) const;

    [[nodiscard]] ProcessIdentifier
    get_neighbour(const ProcessIdentifier &process_identifier, int32_t dl, int32_t di, int32_t dj,
                  bool modulo = false) const;

    [[nodiscard]] IndexRange
    tell_my_range(
            uint32_t process_id_i, uint32_t process_id_j,
            uint32_t n, uint32_t m, int axis
    );

    [[nodiscard]] IndexRange
    tell_my_range(
            uint32_t process_id_l, uint32_t process_id_i, uint32_t process_id_j,
            uint32_t n, uint32_t m, int axis, SplitOverMatrix split_over
    );

    void clear_cache();

private:
    using KeyRange = std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, int, SplitOverMatrix>;
    std::map<KeyRange, IndexRange> cache;
};

std::ostream &operator<<(std::ostream &os, const Topology &metadata);

std::ostream &operator<<(std::ostream &os, const ProcessIdentifier &identifier);

std::string to_string(const ProcessIdentifier &identifier);

bool is(const ProcessIdentifier &id, uint32_t layer, uint32_t row, uint32_t col);
