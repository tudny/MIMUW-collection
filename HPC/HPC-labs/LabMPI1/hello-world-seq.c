/*
 * A template for the 2016 MPI lab at the University of Warsaw.
 * Copyright (C) 2016, Konrad Iwanicki
 * Further modifications by Krzysztof Rzadca 2018
 */
#define _POSIX_TIMERS

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    struct timespec spec;

    clock_gettime(CLOCK_REALTIME, &spec);
    srand(spec.tv_nsec); // use nsec to have a different value across different processes

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int number_of_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    unsigned t = rand() % 5;
    sleep(t);
    printf("Hello world from %d/%d (slept %u s)!\n", rank, number_of_processes, t);

    MPI_Finalize();

    return 0;
}
