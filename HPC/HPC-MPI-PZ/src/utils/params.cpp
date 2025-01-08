#include "params.hpp"
#include <args.hxx>
#include <iostream>
#include "logging.hpp"

#include "fs_utils.hpp"

std::ostream &MatmulParams::operator<<(std::ostream &os, const MatmulParams::AlgorithmType &type) {
    switch (type) {
        case AlgorithmType::SUMMA2D:
            os << "2D";
            break;
        case AlgorithmType::SUMMA3D:
            os << "3D";
            break;
        case AlgorithmType::balanced3D:
            os << "balanced";
            break;
    }
    return os;
}

std::ostream &MatmulParams::operator<<(std::ostream &os, const MatmulParams::Params &params) {
    os << "Algorithm type : " << params.algorithm_type << std::endl;
    os << "Matrix A file  : " << params.sparse_matrix_file_a << std::endl;
    os << "Matrix B file  : " << params.sparse_matrix_file_b << std::endl;
    os << "Print result   : " << params.print_result << std::endl;
    if (params.layers.has_value()) {
        os << "Layers         : " << params.layers.value() << std::endl;
    } else {
        os << "Layers         : not specified" << std::endl;
    }
    if (params.g_value.has_value()) {
        os << "G value        : " << params.g_value.value() << std::endl;
    } else {
        os << "G value        : not specified" << std::endl;
    }
    return os;
}

std::shared_ptr<MatmulParams::Params> MatmulParams::parse_args(int argc, char *argv[]) {

    args::ArgumentParser parser("Matrix multiplication application");
    args::HelpFlag help(
            parser,
            ARG_HELP,
            "Print help message",
            {'h', ARG_HELP}
    );
    args::ValueFlag<std::string> a_file(
            parser,
            ARG_A_FILE,
            "Path to CSR file storing given sparse matrix",
            {'a', ARG_A_FILE},
            args::Options::Required
    );
    args::ValueFlag<std::string> b_file(
            parser,
            ARG_B_FILE,
            "Path to CSR file storing given sparse matrix",
            {'b', ARG_B_FILE},
            args::Options::Required
    );
    args::Flag verbose(
            parser,
            ARG_VERBOSE,
            "Prints the matrix C (the multiplication result) in the row-major order",
            {'v', ARG_VERBOSE}
    );
    args::ValueFlag<std::string> type(
            parser,
            ARG_TYPE,
            "Type specifies which version of the algorithm you have to use. Possible type values are 2D, 3D and balanced.",
            {'t', ARG_TYPE},
            args::Options::Required
    );
    args::ValueFlag<int> layers(
            parser,
            ARG_LAYERS, "Layers specifies the number of layers in the 3D-SUMMA procedure.",
            {'l', ARG_LAYERS}
    );
    args::ValueFlag<double> g_value(
            parser,
            ARG_G_VALUE, "Prints the number of elements in C greater than the g_value.",
            {'g', ARG_G_VALUE}
    );

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help &) {
        std::cout << parser;
        exit(0);
    } catch (args::ParseError &e) {
        std::cerr << "Invalid arguments. Parsing error: ";
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        exit(1);
    } catch (args::ValidationError &e) {
        std::cerr << "Invalid arguments. Validation error: ";
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        exit(1);
    }

    auto params = std::make_shared<Params>();
    params->sparse_matrix_file_a = a_file.Get();
    params->sparse_matrix_file_b = b_file.Get();
    params->print_result = verbose.Get();
    params->algorithm_type = of_string(type.Get());
    if (is_3d(params->algorithm_type)) {
        if (layers) {
            params->layers = layers.Get();
        } else {
            std::cerr << "Layers is required for " << params->algorithm_type << " algorithm";
            exit(1);
        }
    }
    if (g_value) {
        params->g_value = g_value.Get();
    }

    LOGINFO("Parsed parameters");

    return params;
}

MatmulParams::AlgorithmType MatmulParams::of_string(const std::string &str) {
    if (str == "2D") {
        return AlgorithmType::SUMMA2D;
    } else if (str == "3D") {
        return AlgorithmType::SUMMA3D;
    } else if (str == "balanced") {
        return AlgorithmType::balanced3D;
    } else {
        throw std::invalid_argument("Unknown algorithm type: " + str);
    }
}

constexpr bool MatmulParams::is_3d(MatmulParams::AlgorithmType type) {
    return type == AlgorithmType::SUMMA3D || type == AlgorithmType::balanced3D;
}

void MatmulParams::Params::validate_files() const {
    for (const auto &file: {this->sparse_matrix_file_a, this->sparse_matrix_file_b}) {
        if (!file_exists(file)) {
            LOGERROR("File " << file << " does not exist");
            exit(1);
        }
    }
}
