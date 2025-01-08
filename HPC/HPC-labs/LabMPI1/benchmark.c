#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>

#define NUMBER_OF_TESTS 30
#define MAX_SIZE 10000000

double benchmark(int rank, size_t size) {
    double start_time;
    double end_time;

    void *data = malloc(size);

    start_time = MPI_Wtime();

    if (rank == 0) {
        MPI_Send(data, size, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(data, size, MPI_BYTE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(data, size, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(data, size, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }

    end_time = MPI_Wtime();

    return end_time - start_time;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank;
    int number_of_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    assert(number_of_processes == 2);

    for (size_t size = 1; size <= MAX_SIZE; size *= 10) {
        for (int t = 0; t < NUMBER_OF_TESTS; ++t) {
            double time = benchmark(rank, size);
            if (rank == 0) {
                printf("%d %ld %e\n", t, size, time);
            }
        }
    }

    MPI_Finalize();

    return 0;
}
