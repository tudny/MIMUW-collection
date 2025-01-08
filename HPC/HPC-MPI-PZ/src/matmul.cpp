#include "utils/logging.hpp"
#include "utils/params.hpp"
#include "alg/summa2d.hpp"
#include "alg/summa3d.hpp"
#include "alg/summa3d_balanced.hpp"
#include "driver.hpp"
#include "worker.hpp"
#include "utils/meta.hpp"
#include <mpi.h>
#include <algorithm>

SubAB run_setup(const int rank, const int size, const MatmulParams::ParamsPtr &params) {
    if (rank == 0) {
        return run_driver_setup(params, rank, size);
    } else {
        return run_worker_setup(params, rank, size);
    }
}

void run_collect_result_and_print(
        const MatmulParams::ParamsPtr &params,
        const CompressedSparseRowsMatrixPtr &C,
        const int rank,
        const int size
) {
    if (rank == 0) {
        run_driver_collect_result(params, C, rank, size);
    } else {
        run_worker_collect_result(params, C, rank, size);
    }
}

inline CompressedSparseRowsMatrixPtr run_selected_algorithm(
        const MatmulParams::ParamsPtr &params,
        const Topology &topology,
        const CompressedSparseRowsMatrixPtr &A,
        const CompressedSparseRowsMatrixPtr &B,
        int rank, int size
) {
    switch (params->algorithm_type) {
        case MatmulParams::AlgorithmType::SUMMA2D:
            return run_summa2d(params, topology, A, B, rank, size, 0);
        case MatmulParams::AlgorithmType::SUMMA3D:
            return run_summa3d(params, topology, A, B, rank, size);
        case MatmulParams::AlgorithmType::balanced3D:
            return run_summa3d(params, topology, A, B, rank, size);
    }
    __builtin_unreachable();
}

void run_algorithm(
        const MatmulParams::ParamsPtr &params,
        const CompressedSparseRowsMatrixPtr &A,
        const CompressedSparseRowsMatrixPtr &B,
        int rank, int size
) {
    Topology topology{static_cast<uint32_t>(size), params->layers};

    auto C = run_selected_algorithm(params, topology, A, B, rank, size);

    LOGTIME("matmul-algorithm-done-timestamp");

    LOGINFO("Matrix C: " << C->metadata());

    if (params->print_result) {
        run_collect_result_and_print(params, C, rank, size);
    }

    LOGTIME("matmul-g-value-done-timestamp");

    if (auto g = params->g_value) {
        int64_t elements_greater_than_g = std::count_if(
                C->values.begin(), C->values.end(),
                [g](const matrix_type_t &value) { return value > *g; }
        );
        if (*g < 0.0) {
            // if g is negative, we need to count-in all the zeros
            elements_greater_than_g += static_cast<int64_t>(C->n * C->m - C->nnz);
        }
        LOGINFO("Number of elements greater than g: " << elements_greater_than_g);

        int64_t total_elements_greater_than_g;
        MPI_Reduce(&elements_greater_than_g, &total_elements_greater_than_g, 1, MPI_INT64_T, MPI_SUM, 0,
                   MPI_COMM_WORLD);

        if (rank == 0) {
            LOGINFO("Total number of elements greater than g: " << total_elements_greater_than_g);
            std::cout << total_elements_greater_than_g << std::endl;
        }
    }

    LOGTIME("matmul-data-collection-done-timestamp");
}

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    int mpi_rank;
    int mpi_size;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    LOGDEBUG("Starting matmul application on rank " << mpi_rank << " of " << mpi_size);

    try {
        LOGINFO("Starting matmul application as rank " << mpi_rank << " of " << mpi_size);

        LOGTIME("matmul-start-timestamp");
        auto params = MatmulParams::parse_args(argc, argv);
        auto [A, B] = run_setup(mpi_rank, mpi_size, params);
        LOGDEBUG("Matrices A and B are ready for rank " << mpi_rank);

        LOGTIME("matmul-data-distribution-done-timestamp");

        run_algorithm(params, A, B, mpi_rank, mpi_size);

    } catch (const std::exception &e) {
        LOGERROR("Exception caught: " << e.what());
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    LOGINFO("Finished matmul application");

    MPI_Finalize();

    return 0;
}
