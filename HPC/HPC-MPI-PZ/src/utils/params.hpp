#pragma once

#include <iostream>
#include <optional>
#include <memory>

#include "logging.hpp"

namespace MatmulParams {
#define ARG_HELP "help"
#define ARG_A_FILE "a_file"
#define ARG_B_FILE "b_file"
#define ARG_VERBOSE "verbose"
#define ARG_TYPE "type"
#define ARG_LAYERS "layers"
#define ARG_G_VALUE "g_value"

    enum class AlgorithmType {
        SUMMA2D,
        SUMMA3D,
        balanced3D,
    };

    constexpr bool is_3d(AlgorithmType type);

    AlgorithmType of_string(const std::string &str);

    const static constexpr AlgorithmType BASELINE = AlgorithmType::SUMMA2D;

    struct Params {
        std::string sparse_matrix_file_a;
        std::string sparse_matrix_file_b;
        bool print_result = false;
        AlgorithmType algorithm_type = BASELINE;
        std::optional<int> layers = std::nullopt;
        std::optional<double> g_value = std::nullopt;

        void validate_files() const;
    };

    using ParamsPtr = std::shared_ptr<Params>;

    std::ostream &operator<<(std::ostream &os, const AlgorithmType &type);

    std::ostream &operator<<(std::ostream &os, const Params &params);

    ParamsPtr parse_args(int argc, char *argv[]);
}
