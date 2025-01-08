#ifndef GPUGENV_WORKER_H
#define GPUGENV_WORKER_H

#include <vector>
#include "query.h"
#include "../indexer/indexer.h"

void query_index_on_gpu_linear(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
);

void query_index_on_gpu_binary(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
);

void query_index_on_gpu_binary_oneway(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
);

#endif // GPUGENV_WORKER_H
