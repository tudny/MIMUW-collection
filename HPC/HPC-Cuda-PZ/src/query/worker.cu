#include "worker.h"
#include "../logger.h"

#define BLOCK_SIZE 1024

int64_t highest_bit(int64_t x) {
    int64_t log = 0;
    for (; x; x >>= 1, log++);
    return log;
}

__global__ void query_index_on_gpu_linear_kernel(
        const IndexHash *indexes,
        const IndexHash *queries,
        QueryResult *results,
        size_t indexes_size,
        size_t queries_size
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= queries_size) {
        return;
    }

    IndexHash query = queries[idx];
    QueryResult result{NO_RESULT};

    for (size_t i = 0; i < indexes_size; i++) {
        IndexHash index = indexes[i];
        if (index == query) {
            result.offset_idx = i;
            break;
        }
    }

    results[idx] = result;
}

void query_index_on_gpu_linear(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
) {
    void *device_indexes;
    void *device_queries;
    void *device_results;
    size_t indexes_size = indexes->size() * sizeof(IndexHash);
    size_t queries_size = queries.size() * sizeof(IndexHash);
    size_t results_size = queries.size() * sizeof(QueryResult);

    cudaMalloc(&device_indexes, indexes_size);
    cudaMalloc(&device_queries, queries_size);
    cudaMalloc(&device_results, results_size);

    cudaMemcpy(device_indexes, indexes->data(), indexes_size, cudaMemcpyHostToDevice);
    cudaMemcpy(device_queries, queries.data(), queries_size, cudaMemcpyHostToDevice);

    size_t blocks = (queries.size() + BLOCK_SIZE - 1) / BLOCK_SIZE;
    query_index_on_gpu_linear_kernel<<<blocks, BLOCK_SIZE>>>(
            (IndexHash *) device_indexes,
            (IndexHash *) device_queries,
            (QueryResult *) device_results,
            indexes->size(),
            queries.size()
    );

    cudaMemcpy(results.data(), device_results, results_size, cudaMemcpyDeviceToHost);

    cudaFree(device_indexes);
    cudaFree(device_queries);
    cudaFree(device_results);
}

__global__ void query_index_on_gpu_binary_kernel(
        const IndexHash *indexes,
        const IndexHash *queries,
        QueryResult *results,
        size_t indexes_size,
        size_t queries_size
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= queries_size) {
        return;
    }

    IndexHash query = queries[idx];
    QueryResult result = {NO_RESULT};

    int64_t left = -1;
    int64_t right = (int) indexes_size;
    while (left + 1 < right) {
        int64_t mid = (left + right) / 2;
        IndexHash index = indexes[mid];
        if (index < query) {
            left = mid;
        } else if (index > query) {
            right = mid;
        } else {
            result.offset_idx = mid;
            break;
        }
    }

    results[idx] = result;
}

void query_index_on_gpu_binary(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
) {
    void *device_indexes;
    void *device_queries;
    void *device_results;
    size_t indexes_size = indexes->size() * sizeof(IndexHash);
    size_t queries_size = queries.size() * sizeof(IndexHash);
    size_t results_size = queries.size() * sizeof(QueryResult);

    cudaMalloc(&device_indexes, indexes_size);
    cudaMalloc(&device_queries, queries_size);
    cudaMalloc(&device_results, results_size);

    TIMEIT("Index host->GPU", cudaMemcpy(device_indexes, indexes->data(), indexes_size, cudaMemcpyHostToDevice));
    TIMEIT("Query host->GPU", cudaMemcpy(device_queries, queries.data(), queries_size, cudaMemcpyHostToDevice));

    uint64_t blocks = ((uint64_t) queries.size() + BLOCK_SIZE - 1) / BLOCK_SIZE;
    TIMEIT_BEGIN(kernel);
    query_index_on_gpu_binary_kernel<<<blocks, BLOCK_SIZE>>>(
            (IndexHash *) device_indexes,
            (IndexHash *) device_queries,
            (QueryResult *) device_results,
            indexes->size(),
            queries.size()
    );
    TIMEIT_END(kernel);

    TIMEIT("Result GPU->host", cudaMemcpy(results.data(), device_results, results_size, cudaMemcpyDeviceToHost));

    cudaFree(device_indexes);
    cudaFree(device_queries);
    cudaFree(device_results);
}

__global__ void query_index_on_gpu_binary_kernel_oneway(
        const IndexHash *indexes,
        const IndexHash *queries,
        QueryResult *results,
        size_t indexes_size,
        size_t queries_size,
        int64_t log
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= queries_size) {
        return;
    }
    IndexHash query = queries[idx];

    int64_t begin = -1;
    for (int64_t step = 1 << log; step >= 1; step >>= 1) {
        int64_t next = begin + step;
        if (next < (int64_t) indexes_size && (query >= indexes[next])) {
            begin = next;
        }
    }

    results[idx].offset_idx = NO_RESULT;
    if (begin != -1) {
        IndexHash index = indexes[begin];
        if (index == query) {
            results[idx].offset_idx = begin;
        }
    }
}

void query_index_on_gpu_binary_oneway(
        const std::shared_ptr<std::vector<IndexHash>> &indexes,
        const QueryList &queries,
        std::vector<QueryResult> &results
) {
    void *device_indexes;
    void *device_queries;
    void *device_results;
    size_t indexes_size = indexes->size() * sizeof(IndexHash);
    size_t queries_size = queries.size() * sizeof(IndexHash);
    size_t results_size = queries.size() * sizeof(QueryResult);

    cudaMalloc(&device_indexes, indexes_size);
    cudaMalloc(&device_queries, queries_size);
    cudaMalloc(&device_results, results_size);

    TIMEIT("Index host->GPU", cudaMemcpy(device_indexes, indexes->data(), indexes_size, cudaMemcpyHostToDevice));
    TIMEIT("Query host->GPU", cudaMemcpy(device_queries, queries.data(), queries_size, cudaMemcpyHostToDevice));

    int64_t blocks = ((int64_t) queries.size() + BLOCK_SIZE - 1) / BLOCK_SIZE;
    TIMEIT_BEGIN(kernel);
    query_index_on_gpu_binary_kernel_oneway<<<blocks, BLOCK_SIZE>>>(
            (IndexHash *) device_indexes,
            (IndexHash *) device_queries,
            (QueryResult *) device_results,
            indexes->size(),
            queries.size(),
            highest_bit((int) indexes->size())
    );
    TIMEIT_END(kernel);

    TIMEIT("Result GPU->host", cudaMemcpy(results.data(), device_results, results_size, cudaMemcpyDeviceToHost));

    cudaFree(device_indexes);
    cudaFree(device_queries);
    cudaFree(device_results);
}
