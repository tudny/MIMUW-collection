#pragma once

#include "utils/params.hpp"
#include "matrix/matrix.hpp"
#include "utils/meta.hpp"
#include "utils/tags.hpp"

CompressedSparseRowsMatrixPtr recv_sub_matrix_common(
        const Topology &topology,
        [[maybe_unused]] const ProcessIdentifier &my_id,
        ProcessIdentifier from
);

SubAB run_worker_setup(const MatmulParams::ParamsPtr &params, int mpi_rank, int mpi_size);

void run_worker_collect_result(
        const MatmulParams::ParamsPtr &params,
        const CompressedSparseRowsMatrixPtr &C,
        int rank,
        int size
);
