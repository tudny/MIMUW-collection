#include <gtest/gtest.h>
#include <sstream>
#include "../src/matrix/matrix.hpp"
#include "assets.hpp"

constexpr static const char *MATRIX_MOCK_1 =
        "4 5 6 2\n"
        "6.739934628237549 67.19187812303454 -80.69135603689219 2.142363536732276 -4.205497662154855 48.58265581270115 \n"
        "1 4 4 2 4 0 \n"
        "0 2 3 5 6 \n";

constexpr static const char *MATRIX_MOCK_2 =
        "5 6 9 3\n"
        "16.511649451102244 96.08799118298958 -58.76422223629205 75.76079621771837 -75.3909272139254 -46.49289211688652 69.16457404375237 23.33533766460414 -75.50663553830785 \n"
        "0 1 1 4 2 4 1 2 4 \n"
        "0 2 4 4 6 9 \n";

TEST(MatrixSuite, TestMatrix1Load) {
    std::istringstream input(MATRIX_MOCK_1);
    auto matrix = load_matrix(input);

    ASSERT_EQ(matrix->n, 4);
    ASSERT_EQ(matrix->m, 5);
    ASSERT_EQ(matrix->nnz, 6);

    ASSERT_EQ(matrix->values.size(), 6);
    ASSERT_EQ(matrix->col_indices.size(), 6);
    ASSERT_EQ(matrix->row_offsets.size(), 5);

    ASSERT_DOUBLE_EQ(matrix->values[0], 6.739934628237549);
    ASSERT_DOUBLE_EQ(matrix->values[1], 67.19187812303454);
    ASSERT_DOUBLE_EQ(matrix->values[2], -80.69135603689219);
    ASSERT_DOUBLE_EQ(matrix->values[3], 2.142363536732276);
    ASSERT_DOUBLE_EQ(matrix->values[4], -4.205497662154855);
    ASSERT_DOUBLE_EQ(matrix->values[5], 48.58265581270115);

    ASSERT_EQ(matrix->col_indices[0], 1);
    ASSERT_EQ(matrix->col_indices[1], 4);
    ASSERT_EQ(matrix->col_indices[2], 4);
    ASSERT_EQ(matrix->col_indices[3], 2);
    ASSERT_EQ(matrix->col_indices[4], 4);
    ASSERT_EQ(matrix->col_indices[5], 0);

    ASSERT_EQ(matrix->row_offsets[0], 0);
    ASSERT_EQ(matrix->row_offsets[1], 2);
    ASSERT_EQ(matrix->row_offsets[2], 3);
    ASSERT_EQ(matrix->row_offsets[3], 5);
    ASSERT_EQ(matrix->row_offsets[4], 6);
}

TEST(MatrixSuite, TestMatrix2Load) {
    std::istringstream input(MATRIX_MOCK_2);
    auto matrix = load_matrix(input);

    ASSERT_EQ(matrix->n, 5);
    ASSERT_EQ(matrix->m, 6);
    ASSERT_EQ(matrix->nnz, 9);

    ASSERT_EQ(matrix->values.size(), 9);
    ASSERT_EQ(matrix->col_indices.size(), 9);
    ASSERT_EQ(matrix->row_offsets.size(), 6);

    ASSERT_DOUBLE_EQ(matrix->values[0], 16.511649451102244);
    ASSERT_DOUBLE_EQ(matrix->values[1], 96.08799118298958);
    ASSERT_DOUBLE_EQ(matrix->values[2], -58.76422223629205);
    ASSERT_DOUBLE_EQ(matrix->values[3], 75.76079621771837);
    ASSERT_DOUBLE_EQ(matrix->values[4], -75.3909272139254);
    ASSERT_DOUBLE_EQ(matrix->values[5], -46.49289211688652);
    ASSERT_DOUBLE_EQ(matrix->values[6], 69.16457404375237);
    ASSERT_DOUBLE_EQ(matrix->values[7], 23.33533766460414);
    ASSERT_DOUBLE_EQ(matrix->values[8], -75.50663553830785);

    ASSERT_EQ(matrix->col_indices[0], 0);
    ASSERT_EQ(matrix->col_indices[1], 1);
    ASSERT_EQ(matrix->col_indices[2], 1);
    ASSERT_EQ(matrix->col_indices[3], 4);
    ASSERT_EQ(matrix->col_indices[4], 2);
    ASSERT_EQ(matrix->col_indices[5], 4);
    ASSERT_EQ(matrix->col_indices[6], 1);
    ASSERT_EQ(matrix->col_indices[7], 2);
    ASSERT_EQ(matrix->col_indices[8], 4);

    ASSERT_EQ(matrix->row_offsets[0], 0);
    ASSERT_EQ(matrix->row_offsets[1], 2);
    ASSERT_EQ(matrix->row_offsets[2], 4);
    ASSERT_EQ(matrix->row_offsets[3], 4);
    ASSERT_EQ(matrix->row_offsets[4], 6);
    ASSERT_EQ(matrix->row_offsets[5], 9);
}

TEST(MatrixSuite, TestMatrixMultiplyable) {
    std::istringstream input1(MATRIX_MOCK_1);
    auto matrix1 = load_matrix(input1);

    std::istringstream input2(MATRIX_MOCK_2);
    auto matrix2 = load_matrix(input2);

    ASSERT_TRUE(matrix1->multiplyable(*matrix2));
    ASSERT_FALSE(matrix2->multiplyable(*matrix1));
}

TEST(MatrixSuite, TestMatrixShow) {
    std::istringstream input(MATRIX_MOCK_1);
    auto matrix = load_matrix(input);

    std::ostringstream output;
    show(matrix, output);

    auto [n, m, _] = matrix->metadata();
    std::vector<std::vector<matrix_type_t>> expected_result = {
            {0, 6.73993, 0, 0, 67.1919},
            {0, 0, 0, 0, -80.6914},
            {0, 0, 2.14236, 0, -4.2055},
            {48.5827, 0, 0, 0, 0}
    };

    std::cout << output.str();

    std::istringstream real_output(output.str());

    int real_n, real_m;
    real_output >> real_n >> real_m;
    ASSERT_EQ(n, real_n);
    ASSERT_EQ(m, real_m);

    for (auto &row : expected_result) {
        for (auto &expected_value : row) {
            matrix_type_t value;
            real_output >> value;
            ASSERT_NEAR(value, expected_value, 1e-4);
        }
    }
}

TEST(MatrixSuite, TestMatrixSmallMult) {
    auto a_path = assets::fixtures_dir("small-A.csr");
    auto b_path = assets::fixtures_dir("small-B.csr");
    auto expected_c_path = assets::fixtures_dir("small-C.csr");

    auto A = load_matrix(a_path);
    auto B = load_matrix(b_path);
    auto expected_C = load_matrix(expected_c_path);

    CompressedLocalMatrix C;
    local_matmul(*A, *B, C);
    ASSERT_TRUE(compare_crs_and_mapped(expected_C, C));

    auto C_csr = convert_to_csr(C);
    ASSERT_TRUE(compare_crs_and_csr(expected_C, C_csr));
}

TEST(MatrixSuite, TestMatrixMediumMult) {
    auto a_path = assets::fixtures_dir("medium-A.csr");
    auto b_path = assets::fixtures_dir("medium-B.csr");
    auto expected_c_path = assets::fixtures_dir("medium-C.csr");

    auto A = load_matrix(a_path);
    auto B = load_matrix(b_path);
    auto expected_C = load_matrix(expected_c_path);

    CompressedLocalMatrix C;
    local_matmul(*A, *B, C);
    ASSERT_TRUE(compare_crs_and_mapped(expected_C, C));

    auto C_csr = convert_to_csr(C);
    ASSERT_TRUE(compare_crs_and_csr(expected_C, C_csr));
}

TEST(MatrixSuite, TestMatrixBigMult) {
    auto a_path = assets::fixtures_dir("big-A.csr");
    auto b_path = assets::fixtures_dir("big-B.csr");
    auto expected_c_path = assets::fixtures_dir("big-C.csr");

    auto A = load_matrix(a_path);
    auto B = load_matrix(b_path);
    auto expected_C = load_matrix(expected_c_path);

    CompressedLocalMatrix C;
    local_matmul(*A, *B, C);
    ASSERT_TRUE(compare_crs_and_mapped(expected_C, C));

    auto C_csr = convert_to_csr(C);
    ASSERT_TRUE(compare_crs_and_csr(expected_C, C_csr));
}
