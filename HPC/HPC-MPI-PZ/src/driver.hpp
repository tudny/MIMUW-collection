#pragma once

#include "utils/params.hpp"
#include "matrix/matrix.hpp"
#include "utils/tags.hpp"

RequestsPtr send_sub_matrix(
        const Topology &topology,
        const CompressedSparseRowsMatrix &matrix,
        const PackedMatrixMetadata &metadata,
        const ProcessIdentifier &id
);

SubAB run_driver_setup(const MatmulParams::ParamsPtr &params, int mpi_rank, int mpi_size);

void run_driver_collect_result(
        const MatmulParams::ParamsPtr &params,
        const CompressedSparseRowsMatrixPtr &C,
        int rank,
        int size
);
