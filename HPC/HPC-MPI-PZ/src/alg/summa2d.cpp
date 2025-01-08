#include <mpi.h>
#include "summa2d.hpp"
#include "../utils/tags.hpp"

CompressedSparseRowsMatrixPtr broadcast(
        const CompressedSparseRowsMatrixPtr &M,
        uint32_t sending_i,
        uint32_t sending_j,
        int sending_rank_in_new_world,
        int split_key,
        const ProcessIdentifier &id,
        int rank,
        MPI_Comm comm
) {
    const auto &[_, id_i, id_j] = id;
    MPI_Comm partial_comm;

    MPI_Comm_split(comm, split_key, rank, &partial_comm);

    if (sending_i == id_i && sending_j == id_j) {
        uint64_t messages_count_for_data = (M->nnz + SINGLE_MESSAGE_SIZE_ELEMENTS - 1) / SINGLE_MESSAGE_SIZE_ELEMENTS;
        uint64_t messages_count_for_offsets =
                (M->n + 1 + SINGLE_MESSAGE_SIZE_ELEMENTS - 1) / SINGLE_MESSAGE_SIZE_ELEMENTS;

        PackedMatrixMetadata metadata = M->metadata();
        MPI_Bcast(&metadata, sizeof(metadata), MPI_BYTE, sending_rank_in_new_world, partial_comm);

        for (uint64_t i = 0; i < messages_count_for_data; i++) {
            uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
            uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, M->nnz);
            uint64_t size = end - start;
            MPI_Bcast(M->values.data() + start, static_cast<int>(size), MPI_DOUBLE, sending_rank_in_new_world,
                      partial_comm);
        }

        for (uint64_t i = 0; i < messages_count_for_data; i++) {
            uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
            uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, M->nnz);
            uint64_t size = end - start;
            MPI_Bcast(M->col_indices.data() + start, static_cast<int>(size), MPI_UINT32_T, sending_rank_in_new_world,
                      partial_comm);
        }

        for (uint64_t i = 0; i < messages_count_for_offsets; i++) {
            uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
            uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, static_cast<uint64_t>(M->n) + 1L);
            uint64_t size = end - start;
            MPI_Bcast(M->row_offsets.data() + start, static_cast<int>(size), MPI_UINT64_T, sending_rank_in_new_world,
                      partial_comm);
        }

        MPI_Comm_free(&partial_comm);

        return M;
    } else {
        PackedMatrixMetadata metadata{};
        MPI_Bcast(&metadata, sizeof(metadata), MPI_BYTE, sending_rank_in_new_world, partial_comm);

        uint64_t messages_count_for_data =
                (metadata.nnz + SINGLE_MESSAGE_SIZE_ELEMENTS - 1) / SINGLE_MESSAGE_SIZE_ELEMENTS;
        uint64_t messages_count_for_offsets =
                (metadata.n + 1 + SINGLE_MESSAGE_SIZE_ELEMENTS - 1) / SINGLE_MESSAGE_SIZE_ELEMENTS;

        std::vector<matrix_type_t> values(metadata.nnz);
        std::vector<uint32_t> col_indices(metadata.nnz);
        std::vector<uint64_t> row_offsets(metadata.n + 1);

        for (uint64_t i = 0; i < messages_count_for_data; i++) {
            uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
            uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, metadata.nnz);
            uint64_t size = end - start;
            MPI_Bcast(values.data() + start, static_cast<int>(size), MPI_DOUBLE, sending_rank_in_new_world,
                      partial_comm);
        }

        for (uint64_t i = 0; i < messages_count_for_data; i++) {
            uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
            uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, metadata.nnz);
            uint64_t size = end - start;
            MPI_Bcast(col_indices.data() + start, static_cast<int>(size), MPI_UINT32_T, sending_rank_in_new_world,
                      partial_comm);
        }

        for (uint64_t i = 0; i < messages_count_for_offsets; i++) {
            uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
            uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, static_cast<uint64_t>(metadata.n) + 1L);
            uint64_t size = end - start;
            MPI_Bcast(row_offsets.data() + start, static_cast<int>(size), MPI_UINT64_T, sending_rank_in_new_world,
                      partial_comm);
        }

        MPI_Comm_free(&partial_comm);

        // error: cannot bind packed field ‘metadata.PackedMatrixMetadata::n’ to ‘unsigned int&’
        auto n = metadata.n;
        auto m = metadata.m;
        auto nnz = metadata.nnz;
        return std::make_shared<CompressedSparseRowsMatrix>(
                n, m, nnz,
                values, col_indices, row_offsets
        );
    }
}

CompressedLocalMatrixPtr run_summa2d_without_converting(
        const MatmulParams::ParamsPtr &params,
        const Topology &topology,
        const CompressedSparseRowsMatrixPtr &A,
        const CompressedSparseRowsMatrixPtr &B,
        int rank, int size, uint32_t layer
) {
    (void) params;
    (void) size;

    ProcessIdentifier id = topology.get_process_identifier(static_cast<uint32_t>(rank));
    const auto &[id_l, id_i, id_j] = id;

    assert(id_l == layer);

    MPI_Comm layer_comm;
    MPI_Comm_split(MPI_COMM_WORLD, layer, rank, &layer_comm);

    CompressedLocalMatrixPtr C = std::make_shared<CompressedLocalMatrix>();
    for (uint32_t stage = 0; stage < topology.p_c; ++stage) {
        CompressedSparseRowsMatrixPtr A_recv = broadcast(
                A, id_i, stage, static_cast<int>(stage),
                static_cast<int>(id_i), id, rank, layer_comm
        );

        CompressedSparseRowsMatrixPtr B_recv = broadcast(
                B, stage, id_j, static_cast<int>(stage),
                static_cast<int>(id_j), id, rank, layer_comm
        );

        local_matmul(*A_recv, *B_recv, *C);
    }

    return C;
}

CompressedSparseRowsMatrixPtr run_summa2d(
        const MatmulParams::ParamsPtr &params,
        const Topology &topology,
        const CompressedSparseRowsMatrixPtr &A,
        const CompressedSparseRowsMatrixPtr &B,
        int rank, int size, uint32_t layer
) {
    auto C = run_summa2d_without_converting(params, topology, A, B, rank, size, layer);
    return convert_to_csr(*C);
}
