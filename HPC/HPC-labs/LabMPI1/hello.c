#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv); /* intialize the library with parameters caught by the runtime */

    int rank; /* rank of the process */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* get the rank of the process */
    printf("Hello from process %d\n", rank); /* print the rank of the process */

    if (rank == 0) {
        int size; /* number of processes */
        MPI_Comm_size(MPI_COMM_WORLD, &size); /* get the number of processes */
        printf("Number of processes: %d\n", size); /* print the number of processes */
    }

    MPI_Finalize(); /* mark that we've finished communicating */

    return 0;
}
