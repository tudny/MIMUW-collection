#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int number_of_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    if (rank == 0) {
        int number = 1;
        MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&number, 1, MPI_INT, number_of_processes - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Final number: %d\n", number);
    } else {
        int number;
        MPI_Recv(&number, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        number *= rank;
        MPI_Send(&number, 1, MPI_INT, (rank + 1) % number_of_processes, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}
