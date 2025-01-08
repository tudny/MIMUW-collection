#pragma once

#include "../matrix/matrix.hpp"
#include "../utils/params.hpp"
#include "../utils/meta.hpp"

CompressedLocalMatrixPtr run_summa2d_without_converting(
        const MatmulParams::ParamsPtr &params,
        const Topology &topology,
        const CompressedSparseRowsMatrixPtr &A,
        const CompressedSparseRowsMatrixPtr &B,
        int rank, int size, uint32_t layer
);

CompressedSparseRowsMatrixPtr run_summa2d(
        const MatmulParams::ParamsPtr &params,
        const Topology &topology,
        const CompressedSparseRowsMatrixPtr &A,
        const CompressedSparseRowsMatrixPtr &B,
        int rank, int size,
        uint32_t layer
);
