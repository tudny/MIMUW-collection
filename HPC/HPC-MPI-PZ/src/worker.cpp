#include "worker.hpp"
#include "utils/meta.hpp"
#include "utils/tags.hpp"
#include "matrix/matrix.hpp"


CompressedSparseRowsMatrixPtr recv_sub_matrix_common(
        const Topology &topology,
        [[maybe_unused]] const ProcessIdentifier &my_id,
        ProcessIdentifier from
) {
    int from_rank = topology.get_process_rank_int(from);

    LOGDEBUG("Worker " << my_id << " receiving matrix from " << from);

    PackedMatrixMetadata metadata = {};
    MPI_Recv(&metadata, sizeof(metadata),
             MPI_BYTE, from_rank,
             tags::INIT_MATRIX_METADATA, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE
    );

    LOGDEBUG("Worker " << my_id << " received metadata " << metadata);

    CompressedSparseRowsMatrixPtr matrix = std::make_shared<CompressedSparseRowsMatrix>(metadata);

    LOGDEBUG("Worker " << my_id << " created matrix");

    uint64_t messages_count_for_data = (matrix->nnz + SINGLE_MESSAGE_SIZE_ELEMENTS - 1) / SINGLE_MESSAGE_SIZE_ELEMENTS;
    uint64_t messages_count_for_offsets =
            (matrix->n + 1 + SINGLE_MESSAGE_SIZE_ELEMENTS - 1) / SINGLE_MESSAGE_SIZE_ELEMENTS;

    LOGDEBUG("Worker " << my_id << " calculated messages count");

    for (uint64_t i = 0; i < messages_count_for_data; i++) {
        uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
        uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, matrix->nnz);

        MPI_Recv(matrix->values.data() + start, static_cast<int>(end - start),
                 MPI_DOUBLE, from_rank,
                 tags::INIT_MATRIX_VALUES, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE
        );
    }

    LOGDEBUG("Worker " << my_id << " received values");

    for (uint64_t i = 0; i < messages_count_for_data; i++) {
        uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
        uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, matrix->nnz);

        MPI_Recv(matrix->col_indices.data() + start, static_cast<int>(end - start),
                 MPI_UINT32_T, from_rank,
                 tags::INIT_MATRIX_COL_INDICES, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE
        );
    }

    LOGDEBUG("Worker " << my_id << " received col indices");

    for (uint64_t i = 0; i < messages_count_for_offsets; i++) {
        uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
        uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, static_cast<uint64_t>(matrix->n) + 1L);

        MPI_Recv(matrix->row_offsets.data() + start, static_cast<int>(end - start),
                 MPI_UINT64_T, from_rank,
                 tags::INIT_MATRIX_ROW_OFFSETS, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE
        );
    }

    LOGDEBUG("Worker " << my_id << " received row offsets");

    return matrix;
}

CompressedSparseRowsMatrixPtr recv_sub_matrix(
        const Topology &topology,
        [[maybe_unused]] const ProcessIdentifier &id
) {
    assert(id != ProcessIdentifierMASTER);
    return recv_sub_matrix_common(topology, id, ProcessIdentifierMASTER);
}

SubAB run_worker_setup([[maybe_unused]] const MatmulParams::ParamsPtr &params, int mpi_rank, int mpi_size) {
    assert(mpi_rank != 0 && mpi_rank < mpi_size);

    LOGDEBUG("Worker " << mpi_rank << " started");

    Topology topology{static_cast<uint32_t>(mpi_size), params->layers};

    LOGDEBUG("Worker " << mpi_rank << " topology calculated");

    ProcessIdentifier id = topology.get_process_identifier_int(mpi_rank);

    LOGDEBUG("Worker " << mpi_rank << " process identifier calculated");

    CompressedSparseRowsMatrixPtr my_a_matrix = recv_sub_matrix(topology, id);
    LOGDEBUG("Worker " << mpi_rank << " received A matrix");
    CompressedSparseRowsMatrixPtr my_b_matrix = recv_sub_matrix(topology, id);
    LOGDEBUG("Worker " << mpi_rank << " received B matrix");

    return {my_a_matrix, my_b_matrix};
}

void run_worker_collect_result(
        const MatmulParams::ParamsPtr &params,
        const CompressedSparseRowsMatrixPtr &C,
        int mpi_rank,
        int mpi_size
) {
    Topology topology{static_cast<uint32_t>(mpi_size), params->layers};
    ProcessIdentifier id = topology.get_process_identifier_int(mpi_rank);
    (void) id;

    assert(id != ProcessIdentifierMASTER);

    int to_rank = topology.get_process_rank_int(ProcessIdentifierMASTER);

    uint64_t messages_count_for_data = (C->nnz + SINGLE_MESSAGE_SIZE_ELEMENTS - 1) / SINGLE_MESSAGE_SIZE_ELEMENTS;
    uint64_t messages_count_for_offsets =
            (C->n + 1 + SINGLE_MESSAGE_SIZE_ELEMENTS - 1) / SINGLE_MESSAGE_SIZE_ELEMENTS;

    PackedMatrixMetadata metadata = C->metadata();
    MPI_Send(&metadata, sizeof(metadata),
             MPI_BYTE, to_rank,
             tags::INIT_MATRIX_METADATA, MPI_COMM_WORLD
    );

    for (uint64_t i = 0; i < messages_count_for_data; i++) {
        uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
        uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, C->nnz);
        uint64_t size = end - start;
        MPI_Send(C->values.data() + start, static_cast<int>(size),
                 MPI_DOUBLE, to_rank,
                 tags::INIT_MATRIX_VALUES, MPI_COMM_WORLD
        );
    }

    for (uint64_t i = 0; i < messages_count_for_data; i++) {
        uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
        uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, C->nnz);
        uint64_t size = end - start;
        MPI_Send(C->col_indices.data() + start, static_cast<int>(size),
                 MPI_UINT32_T, to_rank,
                 tags::INIT_MATRIX_COL_INDICES, MPI_COMM_WORLD
        );
    }

    for (uint64_t i = 0; i < messages_count_for_offsets; i++) {
        uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
        uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, static_cast<uint64_t>(C->n) + 1L);
        uint64_t size = end - start;
        MPI_Send(C->row_offsets.data() + start, static_cast<int>(size),
                 MPI_UINT64_T, to_rank,
                 tags::INIT_MATRIX_ROW_OFFSETS, MPI_COMM_WORLD
        );
    }
}
