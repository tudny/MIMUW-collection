package com.github.tudny.algo;

import java.util.Scanner;

public class ShortestPath {
    /**
     * Constant value representing infinity in my Floyd Warshall algorithm and Graph representation.
     */
    public static final int INF = 1_000_000_009;

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        Graph graph = readGraph(scanner);
        int[][] shortest = algorithmFloydWarshall(graph);
        answerQuery(shortest, scanner);
    }

    /**
     * Reads graph from scanner.
     * First two numbers should be N and M where N is the number of nodes and M is the number of edges.
     * Then M triplets A, B, L where A is the begin, B is the end and L is the length of the edge.
     * 0 <= A, B < N and 0 <= L < INF = 1_000_000_009
     * @param scanner Scanner for any input providing proper graph.
     * @return Graph representation built according to input.
     */
    public static Graph readGraph(Scanner scanner) {
        int N = scanner.nextInt();
        int M = scanner.nextInt();
        String invalidNodeError = "Number of node should be within range of the graph (%d), but is %d.";
        String invalidLenError = "Length of the edge should be within range (0, INF), but id %d.";

        Graph graph = new Graph(N);

        for (int i = 0; i < M; i++) {
            int A = scanner.nextInt();
            int B = scanner.nextInt();
            int len = scanner.nextInt();

            assert 0 <= A && A < N : String.format(invalidNodeError, N, A);
            assert 0 <= B && B < N : String.format(invalidNodeError, N, B);
            assert 0 <= len && len < INF : String.format(invalidLenError, len);

            graph.addEdge(A, B, len);
        }

        return graph;
    }

    /**
     * Implementation of the Floyd-Warshall algorithm searching shortest paths in graph between all nodes.
     * Complexity: O(N^3) where N is the number of nodes in processed graph.
     * @param graph Graph to be processed by Floyd Warshall algorithm.
     * @return Two-dimensional array storing information about shortest paths between all two nodes in graph.
     */
    public static int[][] algorithmFloydWarshall(Graph graph) {
        int N = graph.getN();
        int[][] shortest = new int[N][N];

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                shortest[i][j] = graph.getEdge(i, j);
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

    /**
     * @param shortest Array with shortest paths between all nodes in graph.
     * @param scanner Scanner for place with queries.
     */
    public static void answerQuery(int[][] shortest, Scanner scanner) {
        int queries = scanner.nextInt();

        while (queries --> 0) {
            int A = scanner.nextInt();
            int B = scanner.nextInt();

            int answer = shortest[A][B];
            System.out.println(answer == INF ? "INF" : answer);
        }
    }

    /**
     * Class representing graph. Edges between nodes are stored in matrix of distances.
     */
    private static class Graph {
        /**
         * Number of nodes in graph.
         */
        private final int N;

        /**
         * Matrix of distance. For each node we store:
         * 0                - for dist[i][i], 0<i<N
         * INF              - if edge doesn't exist
         * value : (0, INF) - if edge exist; value is a weight of the edge.
         * It could have been done using enumerator, but it will certainly slower the algorithm.
         */
        private final int[][] dist;

        /**
         * Creates graph object with N nodes.
         * Initially distances between all nodes are infinity excluding distances from A to A which are 0.
         * @param N Number of nodes in graph.
         */
        public Graph(int N) {
            this.N = N;
            this.dist = new int[N][N];

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    dist[i][j] = INF;
                }
                dist[i][i] = 0;
            }
        }

        /**
         * Adds edge to the graph only if there is no other path between those two nodes or its length is shorter.
         * @param begin Edge beginning
         * @param end Edge ending
         * @param len Edge length
         */
        public void addEdge(int begin, int end, int len) {
            dist[begin][end] = Math.min(dist[begin][end], len);
        }

        /**
         * @param begin Edge beginning
         * @param end Edge ending
         * @return Distance between nodes begin and end or if it doesn't exist then INF const.
         */
        public int getEdge(int begin, int end) {
            return dist[begin][end];
        }

        /**
         * @return Number of nodes in graph.
         */
        public int getN() {
            return N;
        }
    }
}
