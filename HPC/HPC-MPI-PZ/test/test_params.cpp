#include <gtest/gtest.h>

#include "../src/utils/params.hpp"

class Argv {
public:
    Argv(const std::initializer_list<std::string> args) {
        size_t total_buffer_size = 0;
        for (const auto &arg: args) {
            total_buffer_size += arg.size() + 1;
        }

        _buffer = new char[total_buffer_size];
        _argv = new char *[args.size()];
        _size = static_cast<int>(args.size());

        char *current = _buffer;
        size_t i = 0;
        for (const auto &arg: args) {
            _argv[i++] = current;
            std::copy(arg.begin(), arg.end(), current);
            current += arg.size();
            *current = '\0';
            current++;
        }
    }

    ~Argv() {
        delete[] _buffer;
    }

    [[nodiscard]] int size() const { return _size; }

    char **data() { return _argv; }

private:
    char *_buffer;
    char **_argv;
    int _size;
};

TEST(ParamsSuite, MissingArgs) {
    auto argv = Argv({"./matmul"});
    EXPECT_EXIT(
            MatmulParams::parse_args(argv.size(), argv.data()),
            ::testing::ExitedWithCode(1),
            "Invalid arguments"
    );
}

TEST(ParamsSuite, Help) {
    auto argv = Argv({"./matmul", "--help"});
    MatmulParams::parse_args(argv.size(), argv.data());
    EXPECT_EXIT(
            MatmulParams::parse_args(argv.size(), argv.data()),
            ::testing::ExitedWithCode(0),
            ""
    );
}

TEST(ParamsSuite, AllArgsPassed) {
    auto argv = Argv({"./matmul", "-a", "a.txt", "-b", "b.txt", "-v", "-g", "1", "-t", "2D"});
    auto params = MatmulParams::parse_args(argv.size(), argv.data());
    EXPECT_EQ(params->sparse_matrix_file_a, "a.txt");
    EXPECT_EQ(params->sparse_matrix_file_b, "b.txt");
    EXPECT_EQ(params->print_result, true);
    EXPECT_EQ(params->g_value, 1);
    EXPECT_EQ(params->algorithm_type, MatmulParams::AlgorithmType::SUMMA2D);
    ASSERT_FALSE(params->layers.has_value());
}

TEST(ParamsSuite, AllArgsPassed3D) {
    auto argv = Argv({"./matmul", "-a", "a.txt", "-b", "b.txt", "-v", "-g", "1", "-t", "3D", "-l", "2"});
    auto params = MatmulParams::parse_args(argv.size(), argv.data());
    EXPECT_EQ(params->sparse_matrix_file_a, "a.txt");
    EXPECT_EQ(params->sparse_matrix_file_b, "b.txt");
    EXPECT_EQ(params->print_result, true);
    EXPECT_EQ(params->g_value, 1);
    EXPECT_EQ(params->algorithm_type, MatmulParams::AlgorithmType::SUMMA3D);
    ASSERT_TRUE(params->layers.has_value());
}

TEST(ParamsSuite, AllArgsPassedBalanced) {
    auto argv = Argv({"./matmul", "-a", "a.txt", "-b", "b.txt", "-v", "-g", "1", "-t", "balanced", "-l", "2"});
    auto params = MatmulParams::parse_args(argv.size(), argv.data());
    EXPECT_EQ(params->sparse_matrix_file_a, "a.txt");
    EXPECT_EQ(params->sparse_matrix_file_b, "b.txt");
    EXPECT_EQ(params->print_result, true);
    EXPECT_EQ(params->g_value, 1);
    EXPECT_EQ(params->algorithm_type, MatmulParams::AlgorithmType::balanced3D);
    ASSERT_TRUE(params->layers.has_value());
    EXPECT_EQ(params->layers.value(), 2);
}

TEST(ParamsSuite, AllArgsPassedBalancedNoVerbose) {
    auto argv = Argv({"./matmul", "-a", "a.txt", "-b", "b.txt", "-g", "1", "-t", "balanced", "-l", "2"});
    auto params = MatmulParams::parse_args(argv.size(), argv.data());
    EXPECT_EQ(params->sparse_matrix_file_a, "a.txt");
    EXPECT_EQ(params->sparse_matrix_file_b, "b.txt");
    EXPECT_EQ(params->print_result, false);
    EXPECT_EQ(params->g_value, 1);
    EXPECT_EQ(params->algorithm_type, MatmulParams::AlgorithmType::balanced3D);
    ASSERT_TRUE(params->layers.has_value());
    EXPECT_EQ(params->layers.value(), 2);
}

TEST(ParamsSuite, AllArgsPassedBalancedNoVerboseNoG) {
    auto argv = Argv({"./matmul", "-a", "a.txt", "-b", "b.txt", "-t", "balanced", "-l", "2"});
    auto params = MatmulParams::parse_args(argv.size(), argv.data());
    EXPECT_EQ(params->sparse_matrix_file_a, "a.txt");
    EXPECT_EQ(params->sparse_matrix_file_b, "b.txt");
    EXPECT_EQ(params->print_result, false);
    ASSERT_FALSE(params->g_value.has_value());
    EXPECT_EQ(params->algorithm_type, MatmulParams::AlgorithmType::balanced3D);
    ASSERT_TRUE(params->layers.has_value());
    EXPECT_EQ(params->layers.value(), 2);
}

TEST(ParamsSuite, AllArgsPassedBalancedNoVerboseNoGNoL) {
    auto argv = Argv({"./matmul", "-a", "a.txt", "-b", "b.txt", "-t", "balanced"});
    EXPECT_EXIT(MatmulParams::parse_args(argv.size(), argv.data()), ::testing::ExitedWithCode(1),
                "Layers is required for balanced algorithm");
}

TEST(ParamsSuite, AllArgsPassedBalancedNoVerboseNoGNoLNoB) {
    auto argv = Argv({"./matmul", "-a", "a.txt", "-b", "b.txt", "-t", "balanced"});
    EXPECT_EXIT(MatmulParams::parse_args(argv.size(), argv.data()), ::testing::ExitedWithCode(1),
                "Layers is required for balanced algorithm");
}

TEST(ParamsSuite, AllArgsPassed3DNoLayer) {
    auto argv = Argv({"./matmul", "-a", "a.txt", "-b", "b.txt", "-t", "3D"});
    EXPECT_EXIT(MatmulParams::parse_args(argv.size(), argv.data()), ::testing::ExitedWithCode(1),
                "Layers is required for 3D algorithm");
}
