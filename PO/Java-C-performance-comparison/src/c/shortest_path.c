#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define INF 1000000009

struct __Graph;
typedef struct __Graph Graph;

Graph* createGraph(int);
void addEdge(Graph*, int, int, int);
int getEdge(Graph*, int, int);

Graph* readGraph();
int **algorithmFloydWarshall(Graph*);
void answerQuedy(int**);

int main(void) {
    Graph *graph = readGraph();
    int **shortest = algorithmFloydWarshall(graph);
    answerQuedy(shortest);

    return 0;
}

void *safeMalloc(size_t size) {
    void *ptr = malloc(size);

    if (ptr == NULL) {
        fprintf(stderr, "Memory outage.\n");
        exit(1);
    }

    return ptr;
}

int **createTwoDimArray(int X, int Y) {
    int **array = safeMalloc(sizeof(int*) * X);
    for (int i = 0; i < X; i++) {
        array[i] = safeMalloc(sizeof(int) * Y);
    }

    return array;
}

int min(int a, int b) {
    return a < b ? a : b;
}

Graph* readGraph() {
    int N, M, A, B, L;
    if (scanf(" %d %d", &N, &M) != 2)
        exit(1);

    Graph *graph = createGraph(N);

    while (M --> 0) {
        if (scanf(" %d %d %d", &A, &B, &L) != 3)
            exit(1);

        assert(0 <= A && A < N);
        assert(0 <= B && B < N);
        assert(0 <= L && L < INF);

        addEdge(graph, A, B, L);
    }

    return graph;
}

struct __Graph {
    int N; ///< Number of nodes.
    int **dist; ///< Distance matrix.
};

Graph* createGraph(int N) {
    Graph* graph = safeMalloc(sizeof(Graph));

    graph->N = N;
    graph->dist = createTwoDimArray(N, N);
    for (int i = 0; i < N; i++) {

        for (int j = 0; j < N; j++) {
            graph->dist[i][j] = INF;
        }
        graph->dist[i][i] = INF;
    }

    return graph;
}

void addEdge(Graph *graph, int begin, int end, int len) {
    graph->dist[begin][end] = min(graph->dist[begin][end], len);
}

int getEdge(Graph *graph, int begin, int end) {
    return graph->dist[begin][end];
}

int **algorithmFloydWarshall(Graph *graph) {
    int N = graph->N;
    int **shortest = createTwoDimArray(N, N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            shortest[i][j] = getEdge(graph, i, j);
        }
    }

    for (int k = 0; k < N; k++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (shortest[i][j] > shortest[i][k] + shortest[k][j]) {
                    shortest[i][j] = shortest[i][k] + shortest[k][j];
                }
            }
        }
    }

    return shortest;
}

void answerQuedy(int **shortest) {
    int queries, A, B, answer;
    if (scanf(" %d", &queries) != 1)
        exit(1);

    while (queries --> 0) {
        if (scanf(" %d %d", &A, &B) != 2)
            exit(1);
            
        answer = shortest[A][B];

        if (answer == INF)
            printf("INF\n");
        else
            printf("%d\n", answer);
    }
}