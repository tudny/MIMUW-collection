#pragma once

#include "../matrix/matrix.hpp"
#include "../utils/params.hpp"
#include "../utils/meta.hpp"

CompressedSparseRowsMatrixPtr run_summa3d(
        const MatmulParams::ParamsPtr &params,
        const Topology &topology,
        const CompressedSparseRowsMatrixPtr &A,
        const CompressedSparseRowsMatrixPtr &B,
        int rank, int size
);
