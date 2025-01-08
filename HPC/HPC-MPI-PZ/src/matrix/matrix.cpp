#include <fstream>
#include <set>
#include <functional>
#include <iomanip>
#include "matrix.hpp"


CompressedSparseRowsMatrix::CompressedSparseRowsMatrix(
        uint32_t n, uint32_t m, uint64_t nnz,
        const std::vector<matrix_type_t> &values,
        const std::vector<uint32_t> &col_indices,
        const std::vector<uint64_t> &row_offsets
) : n(n), m(m), nnz(nnz), values(values), col_indices(col_indices), row_offsets(row_offsets) {
    assert(values.size() == nnz);
    assert(col_indices.size() == nnz);
    assert(row_offsets.size() == n + 1);

    DO_IN_DEBUG(
            for (uint32_t i = 0l; i < n; i++) {
                assert(row_offsets[i] <= row_offsets[i + 1]);
            }

            for (uint32_t i = 0; i < nnz; i++) {
                assert(col_indices[i] < m);
            }
    );
}

bool PackedMatrixMetadata::multiplyable(const PackedMatrixMetadata &other) const {
    return m == other.n;
}

bool CompressedSparseRowsMatrix::multiplyable(const CompressedSparseRowsMatrix &other) const {
    return m == other.n;
}

PackedMatrixMetadata CompressedSparseRowsMatrix::metadata() const {
    return {n, m, nnz};
}

CompressedSparseRowsMatrix::CompressedSparseRowsMatrix(uint32_t n, uint32_t m, uint64_t nnz) {
    this->n = n;
    this->m = m;
    this->nnz = nnz;
    this->values.resize(nnz);
    this->col_indices.resize(nnz);
    this->row_offsets.resize(n + 1);
}

CompressedSparseRowsMatrix::CompressedSparseRowsMatrix(const PackedMatrixMetadata &metadata) {
    n = metadata.n;
    m = metadata.m;
    nnz = metadata.nnz;
    values.resize(nnz);
    col_indices.resize(nnz);
    row_offsets.resize(n + 1);
}

CompressedSparseRowsMatrix::CompressedSparseRowsMatrix(
        uint32_t n, uint32_t m, uint64_t nnz,
        std::vector<matrix_type_t> &&values,
        std::vector<uint32_t> &&col_indices,
        std::vector<uint64_t> &&row_offsets
) {
    this->n = n;
    this->m = m;
    this->nnz = nnz;
    this->values = std::move(values);
    this->col_indices = std::move(col_indices);
    this->row_offsets = std::move(row_offsets);
}

CompressedSparseRowsMatrixPtr load_matrix(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        LOGERROR("Cannot open file " << filename);
        throw std::runtime_error("Cannot open file " + filename);
    }

    return load_matrix(file);
}

PackedMatrixMetadata load_metadata(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        LOGERROR("Cannot open file " << filename);
        throw std::runtime_error("Cannot open file " + filename);
    }

    uint32_t n, m, nnz;
    file >> n >> m >> nnz;
    return {n, m, nnz};
}

CompressedSparseRowsMatrixPtr load_matrix(std::istream &input_stream) {
    uint32_t n, m;
    uint64_t nnz, max_nnz_in_row;
    input_stream >> n >> m >> nnz >> max_nnz_in_row;

    std::vector<matrix_type_t> values(nnz);
    std::vector<uint32_t> col_indices(nnz);
    std::vector<uint64_t> row_offsets(n + 1);

    for (uint64_t i = 0; i < nnz; i++) {
        input_stream >> values[i];
    }

    for (uint64_t i = 0; i < nnz; i++) {
        input_stream >> col_indices[i];
    }

    for (uint32_t i = 0; i < n + 1; i++) {
        input_stream >> row_offsets[i];
    }

    return std::make_shared<CompressedSparseRowsMatrix>(n, m, nnz, values, col_indices, row_offsets);
}

void local_matmul(const CompressedSparseRowsMatrix &A, const CompressedSparseRowsMatrix &B,
                  CompressedLocalMatrix &result_tuple) {
    assert(A.multiplyable(B));

    auto &[result, n, m] = result_tuple;
    n = A.n;
    m = B.m;

    for (uint32_t offset_idx = 0; offset_idx < A.n; offset_idx++) {
        for (uint64_t row_idx = A.row_offsets[offset_idx]; row_idx < A.row_offsets[offset_idx + 1]; row_idx++) {
            matrix_type_t a_ij = A.values[row_idx];
            uint32_t i = offset_idx;
            uint32_t j = A.col_indices[row_idx];

            for (uint64_t row_idx_b = B.row_offsets[j]; row_idx_b < B.row_offsets[j + 1]; row_idx_b++) {
                matrix_type_t b_jk = B.values[row_idx_b];
                uint32_t k = B.col_indices[row_idx_b];

                result[to_coord(i, k)] += a_ij * b_jk;
            }
        }
    }

    for (auto it = result.cbegin(); it != result.cend();) {
        if (it->second == 0.0) {
            it = result.erase(it);
        } else {
            it++;
        }
    }
}

// Print uncompressed matrix
void show(const CompressedSparseRowsMatrixPtr &matrix, std::ostream &output_stream) {
    output_stream << std::fixed << std::setprecision(6);
    output_stream << matrix->n << " " << matrix->m << std::endl;
    for (uint32_t i = 0; i < matrix->n; i++) {
        for (uint32_t j = 0; j < matrix->m; j++) {
            bool found = false;
            for (uint64_t k = matrix->row_offsets[i]; k < matrix->row_offsets[i + 1]; k++) {
                if (matrix->col_indices[k] == j) {
                    output_stream << matrix->values[k] << " ";
                    found = true;
                    break;
                }
            }
            if (!found) {
                output_stream << "0 ";
            }
        }
        output_stream << std::endl;
    }
}

bool compare_crs_and_mapped(
        const CompressedSparseRowsMatrixPtr &expected,
        const CompressedLocalMatrix &actual,
        matrix_type_t tol
) {
    std::set<CompressedLocalMatrixCoord> visited;
    auto &[actual_map, n, m] = actual;

    for (uint32_t offset_idx = 0; offset_idx < expected->n; offset_idx++) {
        for (uint64_t row_idx = expected->row_offsets[offset_idx];
             row_idx < expected->row_offsets[offset_idx + 1]; row_idx++) {
            matrix_type_t a_ij = expected->values[row_idx];
            uint32_t i = offset_idx;
            uint32_t j = expected->col_indices[row_idx];

            auto it = actual_map.find(to_coord(i, j));
            if (it == actual_map.end()) {
                LOGERROR("Element [" << i << ", " << j << "] not found in actual matrix, expected " << a_ij);
                return false;
            }

            visited.insert(to_coord(i, j));
            if (std::abs(it->second - a_ij) > tol) {
                LOGERROR("Element [" << i << ", " << j << "] has value " << it->second << ", expected " << a_ij);
                LOGERROR("Difference is " << std::abs(it->second - a_ij) << ", tolerance is " << tol);
                return false;
            }
        }
    }

    for (const auto &[coord, value]: actual_map) {
        if (visited.find(coord) == visited.end()) {
            if (std::abs(value) > tol) {
                LOGERROR("Element [" << from_coord_a(coord) << ", " << from_coord_b(coord)
                                     << "] has value "
                                     << value
                                     << ", expected 0");
                return false;
            }
        }
    }

    return true;
}

CompressedSparseRowsMatrixPtr convert_to_csr(const CompressedLocalMatrix &matrix) {
    auto &[actual_map, n, m] = matrix;

    std::vector<matrix_type_t> values;
    std::vector<uint32_t> col_indices;
    std::vector<uint64_t> row_offsets(n + 1);

    std::map<CompressedLocalMatrixCoord, matrix_type_t> actual_map_ordered(actual_map.begin(), actual_map.end());
    for (const auto &[coord, value]: actual_map_ordered) {
        auto i = from_coord_a(coord);
        auto j = from_coord_b(coord);
        values.push_back(value);
        col_indices.push_back(j);
        row_offsets[i + 1]++;
    }

    for (uint32_t i = 1; i < n + 1; i++) {
        row_offsets[i] += row_offsets[i - 1];
    }

    return std::make_shared<CompressedSparseRowsMatrix>(n, m, values.size(), values, col_indices, row_offsets);
}

uint32_t find_idx_in_bounds(const std::vector<uint32_t> &bounds, uint32_t bound) {
    // This can be optimized with binary search
    for (uint32_t i = 0; i < bounds.size() - 1; i++) {
        if (bounds[i] <= bound && bound < bounds[i + 1]) {
            return i;
        }
    }
    throw std::runtime_error("Bound not found, but it should be");
}

std::vector<CompressedSparseRowsMatrixPtr> convert_to_n_csr_s(
        const CompressedLocalMatrix &matrix,
        uint32_t l,
        const std::vector<uint32_t> &bounds
) {
    assert(bounds.size() == l + 1);

    auto &[actual_map, n, _] = matrix;

    std::vector<std::vector<matrix_type_t>> values_s(l);
    std::vector<std::vector<uint32_t>> col_indices_s(l);
    std::vector<std::vector<uint64_t>> row_offsets_s(l, std::vector<uint64_t>(n + 1));

    std::map<CompressedLocalMatrixCoord, matrix_type_t> actual_map_ordered(actual_map.begin(), actual_map.end());
    for (const auto &[coord, value]: actual_map_ordered) {
        auto i = from_coord_a(coord);
        auto j = from_coord_b(coord);
        auto layer = find_idx_in_bounds(bounds, j);
        values_s[layer].push_back(value);
        col_indices_s[layer].push_back(j - bounds[layer]);
        row_offsets_s[layer][i + 1]++;
    }

    std::vector<CompressedSparseRowsMatrixPtr> result;
    result.reserve(l);
    for (uint32_t layer = 0; layer < l; layer++) {
        for (uint32_t i = 1; i < n + 1; i++) {
            row_offsets_s[layer][i] += row_offsets_s[layer][i - 1];
        }
        auto my_m = bounds[layer + 1] - bounds[layer];
        result.push_back(std::make_shared<CompressedSparseRowsMatrix>(
                n, my_m, values_s[layer].size(),
                values_s[layer],
                col_indices_s[layer],
                row_offsets_s[layer]
        ));
    }

    return result;
}

// Assumes that matrices have sorted rows
CompressedSparseRowsMatrixPtr add_sub_csr_s(
        const std::vector<CompressedSparseRowsMatrixPtr> &matrices
) {
    assert(!matrices.empty());
    auto n = matrices[0]->n;
    auto m = matrices[0]->m;

    DO_IN_DEBUG(
            for (const auto &matrix: matrices) {
                assert(matrix->n == n);
                assert(matrix->m == m);
            }
    );

    CompressedSparseRowsMatrixBuilder builder;
    builder.set_n(n).set_m(m);

    std::vector<uint64_t> ptrs(matrices.size());
    std::vector<uint64_t> ptrs_ends(matrices.size());
    for (uint32_t row = 0; row < n; row++) {
        for (uint32_t matrix_idx = 0; matrix_idx < matrices.size(); ++matrix_idx) {
            ptrs[matrix_idx] = matrices[matrix_idx]->row_offsets[row];
            ptrs_ends[matrix_idx] = matrices[matrix_idx]->row_offsets[row + 1];
        }

        while (true) {
            bool all_empty = true;
            for (uint32_t matrix_idx = 0; matrix_idx < matrices.size(); ++matrix_idx) {
                if (ptrs[matrix_idx] < ptrs_ends[matrix_idx]) {
                    all_empty = false;
                    break;
                }
            }
            if (all_empty) {
                break;
            }

            matrix_type_t sum = 0.0;
            for (uint32_t matrix_idx = 0; matrix_idx < matrices.size(); ++matrix_idx) {
                if (ptrs[matrix_idx] < ptrs_ends[matrix_idx]) {
                    sum += matrices[matrix_idx]->values[ptrs[matrix_idx]];
                    ptrs[matrix_idx]++;
                }
            }

            if (std::abs(sum) > EPS) {
                builder.add_value(sum, row);
            }
        }

        builder.new_row();
    }

    return builder.build_ptr();
}

CompressedSparseRowsMatrixPtr add_sub_csr_s_hashing(
        const std::vector<CompressedSparseRowsMatrixPtr> &matrices
) {
    assert(!matrices.empty());

    CompressedLocalMatrix result;
    auto &[result_map, n, m] = result;
    n = matrices[0]->n;
    m = matrices[0]->m;

    for (const auto &matrix: matrices) {
        for (uint32_t i = 0; i < matrix->n; i++) {
            for (uint64_t row_idx = matrix->row_offsets[i]; row_idx < matrix->row_offsets[i + 1]; row_idx++) {
                auto value = matrix->values[row_idx];
                auto j = matrix->col_indices[row_idx];
                result_map[to_coord(i, j)] += value;
            }
        }
    }

    LOGINFO("Restored matrix of size " << n << "x" << m << " with " << result_map.size() << " non-zero elements");

    return convert_to_csr(result);
}

template<typename T>
static bool compare_sub_vectors(
        typename std::vector<T>::const_iterator expected_start,
        typename std::vector<T>::const_iterator expected_end,
        typename std::vector<T>::const_iterator actual_start,
        typename std::vector<T>::const_iterator actual_end,
        std::function<bool(const T &, const T &)> compare
) {
    if (std::distance(expected_start, expected_end) != std::distance(actual_start, actual_end)) {
        return false;
    }

    std::vector<T> expected(expected_start, expected_end);
    std::vector<T> actual(actual_start, actual_end);

    std::sort(expected.begin(), expected.end());
    std::sort(actual.begin(), actual.end());

    return std::equal(expected.begin(), expected.end(), actual.begin(), compare);
}

bool compare_crs_and_csr(const CompressedSparseRowsMatrixPtr &expected, const CompressedSparseRowsMatrixPtr &actual,
                         matrix_type_t tol) {
    if (expected->n != actual->n || expected->m != actual->m || expected->nnz != actual->nnz) {
        LOGERROR("Matrices have different dimensions");
        return false;
    }

    auto n = expected->n;

    for (uint32_t offset_idx = 0; offset_idx < n; offset_idx++) {
        auto expected_row_start = expected->row_offsets[offset_idx];
        auto actual_row_start = actual->row_offsets[offset_idx];
        auto expected_row_end = expected->row_offsets[offset_idx + 1];
        auto actual_row_end = actual->row_offsets[offset_idx + 1];
        if (expected_row_start != actual_row_start || expected_row_end != actual_row_end) {
            LOGERROR("Row " << offset_idx << " has different number of non-zero elements");
            return false;
        }

        if (!compare_sub_vectors<double>(
                expected->values.cbegin() + static_cast<long>(expected_row_start),
                expected->values.cbegin() + static_cast<long>(expected_row_end),
                actual->values.cbegin() + static_cast<long>(actual_row_start),
                actual->values.cbegin() + static_cast<long>(actual_row_end),
                [tol](const double &a, const double &b) { return std::abs(a - b) < tol; }
        )) {
            LOGERROR("Values in row " << offset_idx << " are different");
            return false;
        }

        if (!compare_sub_vectors<uint32_t>(
                expected->col_indices.cbegin() + static_cast<long>(expected_row_start),
                expected->col_indices.cbegin() + static_cast<long>(expected_row_end),
                actual->col_indices.cbegin() + static_cast<long>(actual_row_start),
                actual->col_indices.cbegin() + static_cast<long>(actual_row_end),
                [](const uint32_t &a, const uint32_t &b) { return a == b; }
        )) {
            LOGERROR("Column indices in row " << offset_idx << " are different");
            return false;
        }
    }

    return true;
}

CompressedLocalMatrixCoord to_coord(uint32_t a, uint32_t b) {
    return (static_cast<uint64_t>(a) << sizeof(uint32_t) * 8l) | b;
}

uint32_t from_coord_a(CompressedLocalMatrixCoord coord) {
    return static_cast<uint32_t>(coord >> sizeof(uint32_t) * 8l);
}

uint32_t from_coord_b(CompressedLocalMatrixCoord coord) {
    return static_cast<uint32_t>(coord & ((1l << sizeof(uint32_t) * 8l) - 1l));
}

std::ostream &operator<<(std::ostream &os, const PackedMatrixMetadata &metadata) {
    return os << "PackedMatrixMetadata{n=" << metadata.n << ", m=" << metadata.m << ", nnz=" << metadata.nnz << "}";
}

CompressedSparseRowsMatrixBuilder &CompressedSparseRowsMatrixBuilder::set_n(uint32_t _n) {
    this->n = _n;
    return *this;
}

CompressedSparseRowsMatrixBuilder &CompressedSparseRowsMatrixBuilder::set_m(uint32_t _m) {
    this->m = _m;
    return *this;
}

CompressedSparseRowsMatrixBuilder &
CompressedSparseRowsMatrixBuilder::add_value(matrix_type_t value, uint32_t col_index) {
    values.push_back(value);
    col_indices.push_back(col_index);
    ++nnz;
    return *this;
}

CompressedSparseRowsMatrixBuilder &CompressedSparseRowsMatrixBuilder::new_row() {
    row_offsets.push_back(values.size());
    return *this;
}

CompressedSparseRowsMatrix CompressedSparseRowsMatrixBuilder::build() {
    return {n, m, nnz, values, col_indices, row_offsets};
}

CompressedSparseRowsMatrixPtr CompressedSparseRowsMatrixBuilder::build_ptr() {
    return std::make_shared<CompressedSparseRowsMatrix>(n, m, nnz, values, col_indices, row_offsets);
}

CompressedSparseRowsMatrixBuilder &CompressedSparseRowsMatrixBuilder::extend_n(uint32_t _n) {
    n += _n;
    return *this;
}
