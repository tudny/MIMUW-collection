#include <stdexcept>
#include <cmath>
#include <iostream>
#include <map>
#include "meta.hpp"

ProcessIdentifier ProcessIdentifierMASTER = {0, 0, 0};

std::ostream &operator<<(std::ostream &os, const Topology &metadata) {
    return os << "Topology{p_r=" << metadata.p_r << ", p_c=" << metadata.p_c << ", l=" << metadata.l << "}";
}

std::ostream &operator<<(std::ostream &os, const ProcessIdentifier &identifier) {
    auto &[layer, row, col] = identifier;
    return os << "ProcessIdentifier{layer=" << layer << ", row=" << row << ", col=" << col << "}";
}

bool is(const ProcessIdentifier &id, uint32_t layer, uint32_t row, uint32_t col) {
    auto [l, r, c] = id;
    return l == layer && r == row && c == col;
}

std::string to_string(const ProcessIdentifier &identifier) {
    auto &[layer, row, col] = identifier;
    return "ProcessIdentifier{layer=" + std::to_string(layer) + ", row=" + std::to_string(row) + ", col=" +
           std::to_string(col) + "}";
}

void Topology::validate_matrix_size(uint32_t rows, uint32_t cols) const {
    if (!(rows >= p_r && cols >= p_c)) {
        throw std::runtime_error("Matrix size is not compatible with the topology");
    }
}

ProcessIdentifier Topology::get_process_identifier(uint32_t process_rank) const {
    uint32_t my_layer = process_rank / (p_r * p_c);
    uint32_t my_rank = process_rank % (p_r * p_c);
    uint32_t my_row = my_rank / p_c;
    uint32_t my_col = my_rank % p_c;
    return {my_layer, my_row, my_col};
}

uint32_t Topology::get_process_rank(const ProcessIdentifier &process_identifier) const {
    auto [layer, row, col] = process_identifier;
    return layer * p_r * p_c + row * p_c + col;
}

ProcessIdentifier Topology::get_neighbour(
        const ProcessIdentifier &process_identifier,
        int32_t dl, int32_t di, int32_t dj, bool modulo
) const {
    auto [layer, row, col] = process_identifier;

    int32_t new_layer = static_cast<int32_t>(layer) + dl;
    int32_t new_row = static_cast<int32_t>(row) + di;
    int32_t new_col = static_cast<int32_t>(col) + dj;

    if (!modulo && (new_layer < 0 || new_layer >= static_cast<int32_t>(l) ||
        new_row < 0 || new_row >= static_cast<int32_t>(p_r) ||
        new_col < 0 || new_col >= static_cast<int32_t>(p_c))) {
        throw std::runtime_error("Invalid neighbour");
    }

    if (modulo) {
        new_layer %= static_cast<int32_t>(l);
        new_row %= static_cast<int32_t>(p_r);
        new_col %= static_cast<int32_t>(p_c);
    }

    return {static_cast<uint32_t>(new_layer), static_cast<uint32_t>(new_row), static_cast<uint32_t>(new_col)};
}

Topology::Topology(uint32_t number_of_processes, std::optional<uint32_t> layers) {
    l = layers.value_or(1);
    if (number_of_processes % l != 0) {
        throw std::runtime_error("Number of processes must be divisible by layers");
    }
    uint32_t grid2d = number_of_processes / l;
    auto grid_size = static_cast<uint32_t>(std::sqrt(grid2d));

    if (grid_size * grid_size * l != number_of_processes) {
        throw std::runtime_error("Number of processes must be a square number");
    }

    p_r = grid_size;
    p_c = grid_size;
}

IndexRange Topology::tell_my_range( // NOLINT(*-no-recursion)
        uint32_t process_id_l,
        uint32_t process_id_i,
        uint32_t process_id_j,
        uint32_t n,
        uint32_t m,
        /// 0 - rows, 1 - cols
        int axis,
        SplitOverMatrix split_over
) {
    auto key = std::make_tuple(process_id_l, process_id_i, process_id_j, n, m, axis, split_over);
    auto key_it = cache.find(key);
    if (key_it != cache.end()) {
        return key_it->second;
    }

    if (axis == 0) {
        uint32_t rows_per_process = n / p_r;
        uint32_t remaining_rows = n % p_r;

        uint32_t start = process_id_i * rows_per_process;
        uint32_t end = (process_id_i + 1) * rows_per_process;

        if (process_id_i < remaining_rows) {
            start += process_id_i;
            end += process_id_i + 1;
        } else {
            start += remaining_rows;
            end += remaining_rows;
        }

        // Now we split rows into l layers and take l-th layer
        if (l != 1 and split_over == SplitOverMatrix::B) {
            // i x j part to split into l row chunks
            // We know this range has length `small_n` and we just call same splitter, but with smaller matrix
            auto small_n = end - start;
            auto real_start = start;
            auto [inner_range_start, inner_range_end] = Topology{l * l, std::nullopt}
                    .tell_my_range(process_id_l, 0, small_n, ~0U, 0);
            start = real_start + inner_range_start;
            end = real_start + inner_range_end;
        }

        cache[key] = IndexRange{start, end};
        return IndexRange{start, end};
    } else {
        uint32_t cols_per_process = m / p_c;
        uint32_t remaining_cols = m % p_c;

        uint32_t start = process_id_j * cols_per_process;
        uint32_t end = (process_id_j + 1) * cols_per_process;

        if (process_id_j < remaining_cols) {
            start += process_id_j;
            end += process_id_j + 1;
        } else {
            start += remaining_cols;
            end += remaining_cols;
        }

        // Now we split columns into l layers and take l-th layer
        if (l != 1 and split_over == SplitOverMatrix::A) {
            // ixj part to split into l column chunks
            // We know this range has length `small_m` and we just call same splitter, but with smaller matrix
            auto small_m = end - start;
            auto real_start = start;
            auto [inner_range_start, inner_range_end] = Topology{l * l, std::nullopt}
                    .tell_my_range(0, process_id_l, ~0U, small_m, 1);
            start = real_start + inner_range_start;
            end = real_start + inner_range_end;
        }

        cache[key] = IndexRange{start, end};
        return IndexRange{start, end};
    }
}

IndexRange Topology::tell_my_range( // NOLINT(*-no-recursion)
        uint32_t process_id_i,
        uint32_t process_id_j,
        uint32_t n,
        uint32_t m,
        int axis
) {
    return this->tell_my_range(0, process_id_i, process_id_j, n, m, axis, SplitOverMatrix::I_DONT_CARE);
}

int Topology::get_process_rank_int(const ProcessIdentifier &process_identifier) const {
    return static_cast<int>(get_process_rank(process_identifier));
}

ProcessIdentifier Topology::get_process_identifier_int(int process_rank) const {
    return get_process_identifier(static_cast<uint32_t>(process_rank));
}

void Topology::clear_cache() {
    cache.clear();
}
