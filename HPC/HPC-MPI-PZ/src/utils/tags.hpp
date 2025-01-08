#pragma once

#include <mpi.h>

namespace tags {
    enum TAG {
        INIT_MATRIX_METADATA,
        INIT_MATRIX_VALUES,
        INIT_MATRIX_COL_INDICES,
        INIT_MATRIX_ROW_OFFSETS,
    };
}

using RequestsPtr = std::shared_ptr<std::vector<MPI_Request>>;

/// Single message can contain at most 1M elements
/// For uint64_t, this is 8MB
/// For double, this is 8MB
/// For uint32_t, this is 4MB
constexpr static const uint64_t SINGLE_MESSAGE_SIZE_ELEMENTS = 1024 * 1024;
