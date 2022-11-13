#include <stdio.h>
#include <stdlib.h>

#define MIN_N 0
#define MAX_N 1500

#define MIN_M (MIN_N * MIN_N)
#define MAX_M (MAX_N * 2)

#define MIN_L 0
#define MAX_L 1000000

#define MIN_Q 0
#define MAX_Q 100000

int randomRange(int, int);

int main(int argc, char *argv[]) {

    if (argc != 1) {
        srand(8972); // Generated with Google random
    }
    else {
        int seed = atoi(argv[0]);
        srand(seed);
    }

    int N = randomRange(MIN_N, MAX_N);
    int M = randomRange(MIN_M, MAX_M);
    printf("%d %d\n", N, M);

    while (M --> 0) {
        int A = randomRange(0, N - 1),
            B = randomRange(0, N - 1),
            L = randomRange(MIN_L, MAX_L);

        printf("%d %d %d\n", A, B, L);
    }

    int Q = randomRange(MIN_Q, MAX_Q);
    printf("%d\n", Q);

    while(Q --> 0) {
        int A = randomRange(0, N - 1),
            B = randomRange(0, N - 1);

        printf("%d %d\n", A, B);
    }

    return 0;
}


int randomRange(int a, int b) {
    return a + (rand() % (b - a + 1));
}
