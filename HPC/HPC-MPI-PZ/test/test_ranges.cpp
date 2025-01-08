#include <gtest/gtest.h>
#include <random>
#include "../src/utils/meta.hpp"

void check_ranges(
        uint32_t proc_num_n,
        uint32_t proc_num_m,
        uint32_t n,
        uint32_t m
) {
    Topology topology{proc_num_n * proc_num_m, std::nullopt};

    // Check axis=0
    std::set<uint32_t> lengths;
    for (uint32_t j = 0; j < proc_num_m; j++) {
        auto [start_0, end_0] = topology.tell_my_range(0, j, n, m, 0);
        ASSERT_EQ(start_0, 0);
        lengths.insert(end_0);
        for (uint32_t i = 1; i < proc_num_n; i++) {
            auto [im1_start, im1_end] = topology.tell_my_range(i - 1, j, n, m, 0);
            auto [i_start, i_end] = topology.tell_my_range(i, j, n, m, 0);

            ASSERT_EQ(im1_end, i_start);
            lengths.insert(i_end - i_start);

            // Check if in different axis ranges are the same
            auto [im1_start_1, im1_end_1] = topology.tell_my_range(i - 1, j, n, m, 1);
            auto [i_start_1, i_end_1] = topology.tell_my_range(i, j, n, m, 1);
            ASSERT_EQ(im1_start_1, i_start_1);
            ASSERT_EQ(im1_end_1, i_end_1);
        }
        auto [last_start, last_end] = topology.tell_my_range(proc_num_n - 1, j, n, m, 0);
        ASSERT_EQ(last_end, n);
    }
    ASSERT_EQ(lengths.size(), n % proc_num_n == 0 ? 1 : 2);

    // Check axis=1
    lengths.clear();
    for (uint32_t i = 0; i < proc_num_n; i++) {
        auto [start_1, end_1] = topology.tell_my_range(i, 0, n, m, 1);
        ASSERT_EQ(start_1, 0);
        lengths.insert(end_1);
        for (uint32_t j = 1; j < proc_num_m; j++) {
            auto [jm1_start, jm1_end] = topology.tell_my_range(i, j - 1, n, m, 1);
            auto [j_start, j_end] = topology.tell_my_range(i, j, n, m, 1);

            ASSERT_EQ(jm1_end, j_start);
            lengths.insert(j_end - j_start);

            // Check if in different axis ranges are the same
            auto [jm1_start_0, jm1_end_0] = topology.tell_my_range(i, j - 1, n, m, 0);
            auto [j_start_0, j_end_0] = topology.tell_my_range(i, j, n, m, 0);
            ASSERT_EQ(jm1_start_0, j_start_0);
            ASSERT_EQ(jm1_end_0, j_end_0);
        }
        auto [last_start, last_end] = topology.tell_my_range(i, proc_num_m - 1, n, m, 1);
        ASSERT_EQ(last_end, m);
    }
    ASSERT_EQ(lengths.size(), m % proc_num_m == 0 ? 1 : 2);
}

TEST(RangesSuite, TestRangesRandom) {
    check_ranges(7, 7, 23, 18);

    size_t random_test_count = 100;
    ulong seed = 42;

    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint32_t> proc_num_dist(1, 100);

    for (size_t i = 0; i < random_test_count; i++) {
        uint32_t proc_num_n = proc_num_dist(gen);
        uint32_t proc_num_m = proc_num_n;
        uint32_t n = std::uniform_int_distribution<uint32_t>(proc_num_n, 1000)(gen);
        uint32_t m = std::uniform_int_distribution<uint32_t>(proc_num_m, 1000)(gen);

        check_ranges(proc_num_n, proc_num_m, n, m);
    }
}

void check_layerd_split(
        uint32_t n,
        uint32_t m,
        uint32_t num_proc_n,
        uint32_t num_proc_m,
        uint32_t num_layers
) {
    Topology topology{num_proc_n * num_proc_m * num_layers, num_layers};

    using Key = std::tuple<uint32_t, uint32_t, uint32_t, int, SplitOverMatrix>;
    using Value = IndexRange;
    std::map<Key, Value> cache;

    for (uint32_t proc_i = 0; proc_i < topology.p_r; proc_i++) {
        for (uint32_t proc_j = 0; proc_j < topology.p_c; proc_j++) {
            for (uint32_t proc_l = 0; proc_l < topology.l; proc_l++) {
                auto A_0_range = topology.tell_my_range(proc_l, proc_i, proc_j, n, m, 0, SplitOverMatrix::A);
                auto A_1_range = topology.tell_my_range(proc_l, proc_i, proc_j, n, m, 1, SplitOverMatrix::A);
                auto B_0_range = topology.tell_my_range(proc_l, proc_i, proc_j, n, m, 0, SplitOverMatrix::B);
                auto B_1_range = topology.tell_my_range(proc_l, proc_i, proc_j, n, m, 1, SplitOverMatrix::B);

//                std::cout << "Proc{i=" << proc_i << ", j=" << proc_j << ", l=" << proc_l << "}: "
//                          << "A row: [" << A_0_range.first << ", " << A_0_range.second << "), "
//                          << "A col: [" << A_1_range.first << ", " << A_1_range.second << "), "
//                          << "B row: [" << B_0_range.first << ", " << B_0_range.second << "), "
//                          << "B col: [" << B_1_range.first << ", " << B_1_range.second << ")\n";

                cache[{proc_l, proc_i, proc_j, 0, SplitOverMatrix::A}] = A_0_range;
                cache[{proc_l, proc_i, proc_j, 1, SplitOverMatrix::A}] = A_1_range;
                cache[{proc_l, proc_i, proc_j, 0, SplitOverMatrix::B}] = B_0_range;
                cache[{proc_l, proc_i, proc_j, 1, SplitOverMatrix::B}] = B_1_range;
            }
        }
    }

    for (uint32_t proc_i = 0; proc_i < topology.p_r; proc_i++) {
        std::optional<IndexRange> prev_columns = {};
        std::optional<IndexRange> prev_rows = {};
        for (uint32_t proc_j = 0; proc_j < topology.p_c; proc_j++) {
            std::set<uint64_t> lengths;
            auto [A_all_layers_col_start, A_all_layers_col_end] = topology.tell_my_range(proc_i, proc_j, n, m, 1);
            auto A_all_layers_col_len = A_all_layers_col_end - A_all_layers_col_start;
            for (uint32_t proc_l = 0; proc_l < topology.l; proc_l++) {
                auto [A_col_start, A_col_end] = cache[{proc_l, proc_i, proc_j, 1, SplitOverMatrix::A}];
                auto [A_row_start, A_row_end] = cache[{proc_l, proc_i, proc_j, 0, SplitOverMatrix::A}];
                auto len = A_col_end - A_col_start;
                lengths.insert(len);
                if (prev_columns.has_value()) {
                    ASSERT_EQ(prev_columns->second, A_col_start);
                } else {
                    ASSERT_EQ(A_col_start, 0);
                }
                if (prev_rows.has_value()) {
                    ASSERT_EQ(prev_rows->first, A_row_start);
                    ASSERT_EQ(prev_rows->second, A_row_end);
                }
                prev_columns = {A_col_start, A_col_end};
                prev_rows = {A_row_start, A_row_end};
            }
            ASSERT_EQ(lengths.size(), A_all_layers_col_len % topology.l == 0 ? 1 : 2);
        }
        ASSERT_EQ(prev_columns->second, m);
    }

    {
        std::set<uint64_t> lengths;
        std::optional<IndexRange> prev_rows = {};

        for (uint32_t proc_i = 0; proc_i < topology.p_c; proc_i++) {
            auto [A_row_start, A_row_end] = cache[{0, proc_i, 0, 0, SplitOverMatrix::A}];
            auto len = A_row_end - A_row_start;
            lengths.insert(len);
            if (prev_rows.has_value()) {
                ASSERT_EQ(prev_rows->second, A_row_start);
            } else {
                ASSERT_EQ(A_row_start, 0);
            }
            prev_rows = {A_row_start, A_row_end};
        }
        ASSERT_EQ(prev_rows->second, n);

        ASSERT_EQ(lengths.size(), n % topology.p_r == 0 ? 1 : 2);
    }

    for (uint32_t proc_j = 0; proc_j < topology.p_c; proc_j++) {
        std::optional<IndexRange> prev_columns = {};
        std::optional<IndexRange> prev_rows = {};
        for (uint32_t proc_i = 0; proc_i < topology.p_r; proc_i++) {
            std::set<uint64_t> lengths;
            auto [B_all_layers_row_start, B_all_layers_row_end] = topology.tell_my_range(proc_i, proc_j, n, m, 0);
            auto B_all_layers_row_len = B_all_layers_row_end - B_all_layers_row_start;
            for (uint32_t proc_l = 0; proc_l < topology.l; proc_l++) {
                auto [B_col_start, B_col_end] = cache[{proc_l, proc_i, proc_j, 1, SplitOverMatrix::B}];
                auto [B_row_start, B_row_end] = cache[{proc_l, proc_i, proc_j, 0, SplitOverMatrix::B}];
                auto len = B_row_end - B_row_start;
                lengths.insert(len);
                if (prev_rows.has_value()) {
                    ASSERT_EQ(prev_rows->second, B_row_start);
                } else {
                    ASSERT_EQ(B_row_start, 0);
                }
                if (prev_columns.has_value()) {
                    ASSERT_EQ(prev_columns->first, B_col_start);
                    ASSERT_EQ(prev_columns->second, B_col_end);
                }
                prev_columns = {B_col_start, B_col_end};
                prev_rows = {B_row_start, B_row_end};
            }
            ASSERT_EQ(lengths.size(), B_all_layers_row_len % topology.l == 0 ? 1 : 2);
        }
        ASSERT_EQ(prev_rows->second, n);
    }

    {
        std::set<uint64_t> lengths;
        std::optional<IndexRange> prev_columns = {};

        for (uint32_t proc_j = 0; proc_j < topology.p_r; proc_j++) {
            auto [B_col_start, B_col_end] = cache[{0, 0, proc_j, 1, SplitOverMatrix::B}];
            auto len = B_col_end - B_col_start;
            lengths.insert(len);
            if (prev_columns.has_value()) {
                ASSERT_EQ(prev_columns->second, B_col_start);
            } else {
                ASSERT_EQ(B_col_start, 0);
            }
            prev_columns = {B_col_start, B_col_end};
        }
        ASSERT_EQ(prev_columns->second, m);

        ASSERT_EQ(lengths.size(), m % topology.p_c == 0 ? 1 : 2);
    }
}

TEST(RangesSuite, TestLayerRanges) {
    check_layerd_split(8, 8, 2, 2, 2);
    check_layerd_split(9, 7, 3, 3, 2);
    check_layerd_split(10, 5, 4, 4, 1);
    check_layerd_split(10, 5, 2, 2, 1);

    size_t random_test_count = 100;
    ulong seed = 42;

    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint32_t> proc_num_dist(1, 20);

    for (size_t i = 0; i < random_test_count; i++) {
        uint32_t num_proc_n = proc_num_dist(gen);
        uint32_t num_proc_m = num_proc_n;
        uint32_t n = std::uniform_int_distribution<uint32_t>(num_proc_n, 1000)(gen);
        uint32_t m = std::uniform_int_distribution<uint32_t>(num_proc_m, 1000)(gen);
        uint32_t layers = std::uniform_int_distribution<uint32_t>(0, std::min(n / num_proc_n, m / num_proc_m))(gen) + 1;

        check_layerd_split(n, m, num_proc_n, num_proc_m, layers);
    }
}
