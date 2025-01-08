#pragma once

#include <vector>
#include <map>
#include <unordered_map>
#include <cstdint>
#include "../utils/macros.hpp"
#include "../utils/logging.hpp"
#include <memory>

using matrix_type_t = double;

constexpr static const matrix_type_t EPS = 1e-6;

struct __attribute__((packed)) PackedMatrixMetadata {
    uint32_t n;
    uint32_t m;
    uint64_t nnz;

    [[nodiscard]] bool multiplyable(const PackedMatrixMetadata &other) const;
};

/**
 * -------------
 * | 1 | 0 | 2 |
 * | 0 | 3 | 0 |
 * | 0 | 0 | 4 |
 * | 0 | 5 | 6 |
 * -------------
 * n = 4
 * m = 3
 * nnz = 6
 * values      = {1, 2, 3, 4, 5, 6}
 * col_indices = {0, 2, 1, 2, 1, 2}
 * row_offsets = {0, 2, 3, 4, 6}
 * */
struct CompressedSparseRowsMatrix {
    /// Number of rows
    uint32_t n;
    /// Number of columns
    uint32_t m;
    /// Number of non-zero elements
    uint64_t nnz;
    /// Values of non-zero elements
    std::vector<matrix_type_t> values;
    /// Column indices of non-zero elements
    std::vector<uint32_t> col_indices;
    /// Row offsets
    std::vector<uint64_t> row_offsets;

    CompressedSparseRowsMatrix(
            uint32_t n, uint32_t m, uint64_t nnz,
            const std::vector<matrix_type_t> &values,
            const std::vector<uint32_t> &col_indices,
            const std::vector<uint64_t> &row_offsets
    );

    CompressedSparseRowsMatrix(
            uint32_t n, uint32_t m, uint64_t nnz,
            std::vector<matrix_type_t> &&values,
            std::vector<uint32_t> &&col_indices,
            std::vector<uint64_t> &&row_offsets
    );

    CompressedSparseRowsMatrix(uint32_t n, uint32_t m, uint64_t nnz);

    explicit CompressedSparseRowsMatrix(const PackedMatrixMetadata &metadata);

    [[nodiscard]] bool multiplyable(const CompressedSparseRowsMatrix &other) const;

    [[nodiscard]] PackedMatrixMetadata metadata() const;
};

using CompressedSparseRowsMatrixPtr = std::shared_ptr<CompressedSparseRowsMatrix>;

struct CompressedSparseRowsMatrixBuilder {
private:
    uint32_t n;
    uint32_t m;
    uint64_t nnz = 0;
    std::vector<matrix_type_t> values;
    std::vector<uint32_t> col_indices;
    std::vector<uint64_t> row_offsets = {0};

public:
    CompressedSparseRowsMatrixBuilder &set_n(uint32_t _n);

    CompressedSparseRowsMatrixBuilder &extend_n(uint32_t _n);

    CompressedSparseRowsMatrixBuilder &set_m(uint32_t _m);

    CompressedSparseRowsMatrixBuilder &add_value(matrix_type_t value, uint32_t col_index);

    CompressedSparseRowsMatrixBuilder &new_row();

    CompressedSparseRowsMatrix build();

    CompressedSparseRowsMatrixPtr build_ptr();
};

CompressedSparseRowsMatrixPtr load_matrix(const std::string &filename);

PackedMatrixMetadata load_metadata(const std::string &filename);

CompressedSparseRowsMatrixPtr load_matrix(std::istream &input_stream);

using CompressedLocalMatrixCoord = uint64_t;
using CompressedLocalMatrixMap = std::unordered_map<CompressedLocalMatrixCoord, matrix_type_t>;
using CompressedLocalMatrix = std::tuple<CompressedLocalMatrixMap, uint32_t, uint32_t>;
using CompressedLocalMatrixPtr = std::shared_ptr<CompressedLocalMatrix>;

CompressedLocalMatrixCoord to_coord(uint32_t a, uint32_t b);

uint32_t from_coord_a(CompressedLocalMatrixCoord coord);

uint32_t from_coord_b(CompressedLocalMatrixCoord coord);

void
local_matmul(const CompressedSparseRowsMatrix &A, const CompressedSparseRowsMatrix &B, CompressedLocalMatrix &result);

void show(const CompressedSparseRowsMatrixPtr &matrix, std::ostream &output_stream = std::cout);

bool compare_crs_and_mapped(
        const CompressedSparseRowsMatrixPtr &expected,
        const CompressedLocalMatrix &actual,
        matrix_type_t tol = EPS
);

bool compare_crs_and_csr(
        const CompressedSparseRowsMatrixPtr &expected,
        const CompressedSparseRowsMatrixPtr &actual,
        matrix_type_t tol = EPS
);

CompressedSparseRowsMatrixPtr convert_to_csr(const CompressedLocalMatrix &matrix);


std::vector<CompressedSparseRowsMatrixPtr> convert_to_n_csr_s(
        const CompressedLocalMatrix &matrix,
        uint32_t l,
        const std::vector<uint32_t> &bounds
);

CompressedSparseRowsMatrixPtr add_sub_csr_s(
        const std::vector<CompressedSparseRowsMatrixPtr> &matrices
);

CompressedSparseRowsMatrixPtr add_sub_csr_s_hashing(
        const std::vector<CompressedSparseRowsMatrixPtr> &matrices
);

std::ostream &operator<<(std::ostream &os, const PackedMatrixMetadata &metadata);

using SubAB = std::pair<CompressedSparseRowsMatrixPtr, CompressedSparseRowsMatrixPtr>;
