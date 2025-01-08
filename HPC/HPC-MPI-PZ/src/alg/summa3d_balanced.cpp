#include "summa3d.hpp"

CompressedSparseRowsMatrixPtr run_summa3d_balanced(
        const MatmulParams::ParamsPtr &params,
        const Topology &topology,
        const CompressedSparseRowsMatrixPtr &A,
        const CompressedSparseRowsMatrixPtr &B,
        int rank, int size
) {
    (void) params;
    (void) topology;
    (void) A;
    (void) B;
    (void) rank;
    (void) size;

    throw std::runtime_error("Not implemented yet");
}
