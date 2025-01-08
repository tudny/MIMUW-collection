#include <mpi.h>
#include <csignal>
#include "matrix/matrix.hpp"
#include "utils/meta.hpp"
#include "utils/logging.hpp"
#include "utils/tags.hpp"
#include "driver.hpp"
#include "worker.hpp"

RequestsPtr send_sub_matrix(
        const Topology &topology,
        const CompressedSparseRowsMatrix &matrix,
        const PackedMatrixMetadata &metadata,
        const ProcessIdentifier &id
) {
    int go_rank = topology.get_process_rank_int(id);

    uint64_t messages_count_for_data = (matrix.nnz + SINGLE_MESSAGE_SIZE_ELEMENTS - 1) / SINGLE_MESSAGE_SIZE_ELEMENTS;
    uint64_t messages_count_for_offsets =
            (matrix.n + 1 + SINGLE_MESSAGE_SIZE_ELEMENTS - 1) / SINGLE_MESSAGE_SIZE_ELEMENTS;
    // One message for metadata, one for values, one for col_indices, one for row_offsets
    uint64_t messages_count_total = 1 + 2 * messages_count_for_data + messages_count_for_offsets;

    RequestsPtr requests = std::make_shared<std::vector<MPI_Request>>(messages_count_total);
    uint64_t send_messages_count = 0;

    MPI_Isend(&metadata, sizeof(metadata),
              MPI_BYTE, go_rank,
              tags::INIT_MATRIX_METADATA, MPI_COMM_WORLD,
              &(*requests)[send_messages_count++]
    );

    for (uint64_t i = 0; i < messages_count_for_data; i++) {
        uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
        uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, matrix.nnz);
        uint64_t size = end - start;
        MPI_Isend(matrix.values.data() + start, static_cast<int>(size),
                  MPI_DOUBLE, go_rank,
                  tags::INIT_MATRIX_VALUES, MPI_COMM_WORLD,
                  &(*requests)[send_messages_count++]
        );
    }

    for (uint64_t i = 0; i < messages_count_for_data; i++) {
        uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
        uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, matrix.nnz);
        uint64_t size = end - start;
        MPI_Isend(matrix.col_indices.data() + start, static_cast<int>(size),
                  MPI_UINT32_T, go_rank,
                  tags::INIT_MATRIX_COL_INDICES, MPI_COMM_WORLD,
                  &(*requests)[send_messages_count++]
        );
    }

    for (uint64_t i = 0; i < messages_count_for_offsets; i++) {
        uint64_t start = i * SINGLE_MESSAGE_SIZE_ELEMENTS;
        uint64_t end = std::min((i + 1) * SINGLE_MESSAGE_SIZE_ELEMENTS, static_cast<uint64_t>(matrix.n) + 1L);
        uint64_t size = end - start;
        MPI_Isend(matrix.row_offsets.data() + start, static_cast<int>(size),
                  MPI_UINT64_T, go_rank,
                  tags::INIT_MATRIX_ROW_OFFSETS, MPI_COMM_WORLD,
                  &(*requests)[send_messages_count++]
        );
    }

    assert(send_messages_count == messages_count_total);

    return requests;
}

uint32_t find_idx(const std::vector<uint32_t> &sorted_values, uint32_t value) {
    // TODO: can be optimized
    for (uint32_t i = 0; i < sorted_values.size(); i++) {
        if (sorted_values[i] > value) {
            return i;
        }
    }
    throw std::runtime_error("Value not found");
}

CompressedSparseRowsMatrixPtr split_matrix_and_send(
        Topology &topology,
        const CompressedSparseRowsMatrix &matrix,
        SplitOverMatrix split_matrix
) {
    CompressedSparseRowsMatrixPtr result;

    uint32_t real_proc_over_r = topology.p_r;
    uint32_t real_proc_over_c = topology.p_c;
    switch (split_matrix) {
        case A:
            real_proc_over_c *= topology.l;
            break;
        case B:
            real_proc_over_r *= topology.l;
            break;
        case I_DONT_CARE:
            break;
    }

    LOGINFO("Splitting matrix " << split_matrix << " over grid " << real_proc_over_r << "x" << real_proc_over_c);

    uint32_t current_layer = 0;
    uint32_t real_proc_row = 0;
    uint32_t real_proc_column = 0;
    for (uint32_t process_row = 0; process_row < real_proc_over_r; ++process_row) {
        LOGDEBUG("Splitting row " << process_row);
        if (split_matrix == A) {
            real_proc_row = process_row;
        }
        LOGDEBUG("Calling tell_my_range for layer " << current_layer << " row " << real_proc_row << " column " << 0);
        auto [row_range_start, row_range_end] = topology.tell_my_range(current_layer, real_proc_row, 0, matrix.n,
                                                                       matrix.m, 0, split_matrix);
        LOGDEBUG("Row range: " << row_range_start << " " << row_range_end);
        std::vector<CompressedSparseRowsMatrixBuilder> sub_matrices(real_proc_over_c);
        std::vector<uint32_t> column_ranges_begins(real_proc_over_c + 1);
        LOGDEBUG("Created sub-matrices: " << sub_matrices.size() << " and column ranges: "
                                          << column_ranges_begins.size());
        if (split_matrix == A) {
            current_layer = 0;
        }
        real_proc_column = 0;
        for (uint32_t process_column = 0; process_column < real_proc_over_c; ++process_column) {
            LOGDEBUG("Matrix shape " << matrix.n << "x" << matrix.m);
            LOGDEBUG("Calling tell_my_range for layer " << current_layer << " row " << 0 << " column "
                                                        << real_proc_column);
            auto [column_range_start, column_range_end] =
                    topology.tell_my_range(current_layer, 0, real_proc_column, matrix.n, matrix.m, 1, split_matrix);
            LOGDEBUG("Column range: " << column_range_start << " " << column_range_end);
            column_ranges_begins[process_column] = column_range_start;
            sub_matrices[process_column].set_n(row_range_end - row_range_start);
            sub_matrices[process_column].set_m(column_range_end - column_range_start);
            if (split_matrix == A) {
                current_layer = (current_layer + 1) % topology.l;
                if (current_layer == 0) {
                    real_proc_column++;
                }
            } else {
                ++real_proc_column;
            }
        }
        column_ranges_begins[real_proc_over_c] = matrix.m;

        for (uint32_t offset_idx = row_range_start; offset_idx < row_range_end; offset_idx++) {
            for (uint64_t row_idx = matrix.row_offsets[offset_idx];
                 row_idx < matrix.row_offsets[offset_idx + 1]; row_idx++) {
                matrix_type_t a_ij = matrix.values[row_idx];
                uint32_t j = matrix.col_indices[row_idx];
                auto process_column = find_idx(column_ranges_begins, j) - 1;
                uint32_t local_j = j - column_ranges_begins[process_column];
                sub_matrices[process_column].add_value(a_ij, local_j);
            }
            for (auto &sub_matrix: sub_matrices) {
                sub_matrix.new_row();
            }
        }

        std::vector<RequestsPtr> all_requests;
        std::vector<CompressedSparseRowsMatrix> sub_matrices_built;
        std::vector<PackedMatrixMetadata> sub_matrices_metadata;
        sub_matrices_built.reserve(real_proc_over_c);
        sub_matrices_metadata.reserve(real_proc_over_c);

        // TODO: fix layer
        uint32_t process_column = 0;
        if (split_matrix == A) {
            current_layer = 0;
        }
        for (auto &sub_matrix_builder: sub_matrices) {
            sub_matrices_built.push_back(sub_matrix_builder.build());
            sub_matrices_metadata.push_back(sub_matrices_built.back().metadata());
            auto &sub_matrix = sub_matrices_built.back();
            auto &sub_matrix_metadata = sub_matrices_metadata.back();
            ProcessIdentifier who_to_send = {current_layer, real_proc_row, process_column};

            LOGDEBUG("I'm sending to " << who_to_send << " and their rank is "
                                       << topology.get_process_rank_int(who_to_send));

            if (who_to_send == ProcessIdentifierMASTER) {
                result = std::make_shared<CompressedSparseRowsMatrix>(sub_matrix);
            } else {
                auto requests = send_sub_matrix(topology, sub_matrix, sub_matrix_metadata, who_to_send);
                all_requests.push_back(requests);
            }

            if (split_matrix == A) {
                current_layer = (current_layer + 1) % topology.l;
                if (current_layer == 0) {
                    process_column++;
                }
            } else {
                process_column++;
            }
        }

        for (auto &requests: all_requests) {
            MPI_Waitall(static_cast<int>(requests->size()), requests->data(), MPI_STATUSES_IGNORE);
        }

        // Current layer is incremented only if B is split, as layers make more sub-rows
        if (split_matrix == B) {
            current_layer = (current_layer + 1) % topology.l;
            if (current_layer == 0) {
                real_proc_row++;
            }
        }
    }

    return result;
}

SubAB run_driver_setup(const MatmulParams::ParamsPtr &params, int mpi_rank, int mpi_size) {
    assert(mpi_rank == 0);

    LOGDEBUG("Driver setup");

    params->validate_files();
    LOGINFO(*params << std::endl);

    auto A_meta = load_metadata(params->sparse_matrix_file_a);
    auto B_meta = load_metadata(params->sparse_matrix_file_b);
    if (!A_meta.multiplyable(B_meta)) {
        LOGERROR("Matrices are not multiplyable");
        throw std::runtime_error("Matrices are not multiplyable");
    }

    LOGDEBUG("Matrices are multiplyable");

    Topology topology{static_cast<uint32_t>(mpi_size), params->layers};
    topology.validate_matrix_size(A_meta.n, A_meta.m);
    topology.validate_matrix_size(B_meta.n, B_meta.m);
    LOGINFO("Deduced topology: " << topology);

    ProcessIdentifier id = topology.get_process_identifier_int(mpi_rank);
    (void) id;
    assert(id == ProcessIdentifierMASTER);

    auto A = load_matrix(params->sparse_matrix_file_a);
    LOGDEBUG("Matrix A loaded");
    auto my_a_matrix = split_matrix_and_send(topology, *A, SplitOverMatrix::A);
    LOGDEBUG("Matrix A split and sent");
    A.reset();
    LOGDEBUG("Matrix A released");

    auto B = load_matrix(params->sparse_matrix_file_b);
    LOGDEBUG("Matrix B loaded");
    auto my_b_matrix = split_matrix_and_send(topology, *B, SplitOverMatrix::B);
    LOGDEBUG("Matrix B split and sent");
    B.reset();
    LOGDEBUG("Matrix B released");

    return {my_a_matrix, my_b_matrix};
}

void extend_builder_with(
        CompressedSparseRowsMatrixBuilder &builder,
        std::vector<CompressedSparseRowsMatrixPtr> &matrices
) {
    DO_IN_DEBUG(
            for (auto &matrix: matrices) {
                assert(matrix->n == matrices[0]->n);
            }
    );

    assert(!matrices.empty());
    auto common_sub_n = matrices[0]->n;
    uint32_t whole_m = 0UL;
    for (auto &matrix: matrices) {
        whole_m += matrix->m;
    }
    builder.extend_n(common_sub_n);
    builder.set_m(whole_m);

    for (uint32_t current_row = 0; current_row < common_sub_n; current_row++) {
        uint32_t current_column_offset = 0;
        for (auto &matrix: matrices) {
            for (uint64_t row_idx = matrix->row_offsets[current_row];
                 row_idx < matrix->row_offsets[current_row + 1]; row_idx++) {
                matrix_type_t a_ij = matrix->values[row_idx];
                uint32_t j = matrix->col_indices[row_idx] + current_column_offset;
                builder.add_value(a_ij, j);
            }
            current_column_offset += matrix->m;
        }
        builder.new_row();
    }
}

void run_driver_collect_result(
        const MatmulParams::ParamsPtr &params,
        const CompressedSparseRowsMatrixPtr &my_C,
        [[maybe_unused]] const int mpi_rank,
        const int mpi_size
) {
    assert(mpi_rank == 0);

    CompressedSparseRowsMatrixBuilder C_builder;
    C_builder.set_n(0);

    Topology topology{static_cast<uint32_t>(mpi_size), params->layers};

    uint64_t real_column_width = topology.p_c * topology.l;

    for (uint32_t process_row = 0; process_row < topology.p_r; ++process_row) {
        uint32_t current_layer = 0;
        uint32_t current_column = 0;
        std::vector<CompressedSparseRowsMatrixPtr> sub_matrices(real_column_width);
        for (uint32_t process_column = 0; process_column < real_column_width; ++process_column) {
            ProcessIdentifier who_sends = {current_layer, process_row, current_column};
            ++current_layer;
            if (current_layer == topology.l) {
                current_layer = 0;
                ++current_column;
            }
            if (who_sends == ProcessIdentifierMASTER) {
                sub_matrices[process_column] = my_C;
                continue;
            }
            auto C = recv_sub_matrix_common(topology, ProcessIdentifierMASTER, who_sends);
            LOGINFO("Received matrix C from " << who_sends << " with shape " << C->metadata());
            sub_matrices[process_column] = C;
        }
        extend_builder_with(C_builder, sub_matrices);
    }

    auto C = C_builder.build_ptr();
    LOGINFO("Matrix C: " << C->metadata());
    show(C, std::cout);
}
