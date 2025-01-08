#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#ifndef RADIUS
#define RADIUS -1
#endif

#ifndef NUM_ELEMENTS
#define NUM_ELEMENTS -1
#endif

#define MAX_NUM 2500

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

static void handleError(cudaError_t err, const char *file, int line)
{
    if (err != cudaSuccess)
    {
        printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
        exit(EXIT_FAILURE);
    }
}
#define cudaCheck(err) (handleError(err, __FILE__, __LINE__))

__global__ void stencil_1d(int *in, int *out)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i > NUM_ELEMENTS)
        return;
    int start = MAX(i - RADIUS, 0);
    int end = MIN(i + RADIUS, NUM_ELEMENTS - 1);
    int result = 0;
    for (int k = start; k <= end; ++k)
    {
        result += in[k];
    }
    out[i] = result;
}

void cpu_slow_stencil_1d(int *in, int *out)
{
    #pragma omp parallel for
    for (int i = 0; i < NUM_ELEMENTS; ++i)
    {
        out[i] = 0;
        int start = MAX(i - RADIUS, 0);
        int end = MIN(i + RADIUS, NUM_ELEMENTS - 1);
        for (int k = start; k <= end; ++k)
        {
            out[i] += in[k];
        }
    }
}

void cpu_stencil_1d(int *in, int *out)
{
    int current_sum = 0;
    for (size_t i = 0; i < RADIUS; ++i)
    {
        current_sum += in[i];
    }

    int idxToRemove = -1 - RADIUS - 1;
    int idxToAdd = -1 + RADIUS;
    for (size_t i = 0; i < NUM_ELEMENTS; ++i)
    {
        ++idxToRemove;
        ++idxToAdd;
        if (idxToRemove >= 0)
        {
            current_sum -= in[idxToRemove];
        }
        if (idxToAdd < NUM_ELEMENTS)
        {
            current_sum += in[idxToAdd];
        }
        out[i] = current_sum;
    }
}

int main()
{
    // ************
    // MEMORY SETUP
    // ************

    int *in = (int *)malloc(NUM_ELEMENTS * sizeof(int));
    int *out = (int *)malloc(NUM_ELEMENTS * sizeof(int));
    int *out_cuda_but_on_host = (int *)malloc(NUM_ELEMENTS * sizeof(int));

    for (size_t i = 0; i < NUM_ELEMENTS; ++i)
    {
        in[i] = rand() % MAX_NUM;
    }

    int *cuda_in;
    int *cuda_out;

    // ************
    // SETUP CUDA MEMORY
    // ************

    cudaEvent_t start, stop, start_copy, end_copy, start_kernel, end_kernel, start_restore, end_restore;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventCreate(&start_copy);
    cudaEventCreate(&end_copy);
    cudaEventCreate(&start_kernel);
    cudaEventCreate(&end_kernel);
    cudaEventCreate(&start_restore);
    cudaEventCreate(&end_restore);

    cudaEventRecord(start, 0);

    cudaEventRecord(start_copy, 0);
    cudaCheck(cudaMalloc((void **)&cuda_in, NUM_ELEMENTS * sizeof(int)));
    cudaCheck(cudaMalloc((void **)&cuda_out, NUM_ELEMENTS * sizeof(int)));

    // ************
    // COPY TO GPU
    // ************

    cudaCheck(cudaMemcpy(cuda_in,
                         in,
                         NUM_ELEMENTS * sizeof(int),
                         cudaMemcpyHostToDevice));

    cudaEventRecord(end_copy, 0);

    // ************
    // EXECUTE KERNEL
    // ************

    int blockSize = 256;
    int gridSize = (NUM_ELEMENTS + blockSize - 1) / blockSize;

    cudaEventRecord(start_kernel, 0);
    stencil_1d<<<gridSize, blockSize>>>(cuda_in, cuda_out);
    cudaEventRecord(end_kernel, 0);

    cudaCheck(cudaPeekAtLastError());

    // ************
    // COPY RESULT BACK
    // ************

    cudaEventRecord(start_restore, 0);
    cudaCheck(cudaMemcpy(out_cuda_but_on_host,
                         cuda_out,
                         NUM_ELEMENTS * sizeof(int),
                         cudaMemcpyDeviceToHost));
    cudaEventRecord(end_restore, 0);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float elapsedTime;
    float elapsedCopy;
    float elapsedKernel;
    float elapsedRestore;
    cudaEventElapsedTime(&elapsedTime, start, stop);
    cudaEventElapsedTime(&elapsedCopy, start_copy, end_copy);
    cudaEventElapsedTime(&elapsedKernel, start_kernel, end_kernel);
    cudaEventElapsedTime(&elapsedRestore, start_restore, end_restore);
    printf("Total GPU execution time:  %3.1f ms\n", elapsedTime);
    printf("Total GPU init copy time:  %3.1f ms\n", elapsedCopy);
    printf("Total GPU kernel exe time: %3.1f ms\n", elapsedKernel);
    printf("Total GPU restore time:    %3.1f ms\n", elapsedRestore);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    cudaEventDestroy(start_copy);
    cudaEventDestroy(end_copy);
    cudaEventDestroy(start_kernel);
    cudaEventDestroy(end_kernel);
    cudaEventDestroy(start_restore);
    cudaEventDestroy(end_restore);

    // ************
    // FREE CUDA MEMORY
    // ************

    cudaCheck(cudaFree(cuda_in));
    cudaCheck(cudaFree(cuda_out));

    // ************
    // RUN CPU
    // ************

    struct timespec cpu_start, cpu_stop;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_start);

    cpu_slow_stencil_1d(in, out);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_stop);
    double result = (cpu_stop.tv_sec - cpu_start.tv_sec) * 1e3 + (cpu_stop.tv_nsec - cpu_start.tv_nsec) / 1e6;
    printf("CPU execution time:  %3.1f ms\n", result);


    // ************
    // RUN CPU FAST
    // ************

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_start);

    cpu_stencil_1d(in, out);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_stop);
    result = (cpu_stop.tv_sec - cpu_start.tv_sec) * 1e3 + (cpu_stop.tv_nsec - cpu_start.tv_nsec) / 1e6;
    printf("CPU FAST execution time:  %3.1f ms\n", result);


    // ************
    // SANITY CHECK
    // ************

    for (size_t i = 0; i < NUM_ELEMENTS; ++i)
    {
        if (out[i] != out_cuda_but_on_host[i])
        {
            printf("FAIL: %ld\n", i);
            printf("Expected %d, but got %d\n", out[i], out_cuda_but_on_host[i]);
            break;
        }
    }

    free(in);
    free(out);
    free(out_cuda_but_on_host);

    return 0;
}
