#ifndef GPUGENV_QUERY_H
#define GPUGENV_QUERY_H

#include <vector>
#include "../params.h"
#include "../indexer/indexer.h"

#define NO_RESULT (~0ULL)

struct __attribute__((packed)) QueryResult {
    uint64_t offset_idx = NO_RESULT;
};

using QueryList = std::vector<IndexHash>;

using query_resolution = void (*)(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
);

void query(const QueryParams &params);

#endif // GPUGENV_QUERY_H
