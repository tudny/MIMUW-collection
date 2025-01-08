/*
 * A template for the 2019 MPI lab at the University of Warsaw.
 * Copyright (C) 2016, Konrad Iwanicki.
 * Refactoring 2019, Łukasz Rączkowski
 */

#include <cassert>
#include <mpi.h>
#include "graph-base.h"
#include "graph-utils.h"

int getFirstGraphRowOfProcess(int numVertices, int numProcesses, int myRank) {

    int rowPerProcess = numVertices / numProcesses;
    int rowsLeft = numVertices - rowPerProcess * numProcesses;

    int begin = -1;

    if (myRank < rowsLeft) {
        begin = myRank * (rowPerProcess + 1);
    } else {
        begin = rowsLeft * (rowPerProcess + 1) + (myRank - rowsLeft) * rowPerProcess;
    }

    return begin;
}

Graph* createAndDistributeGraph(int numVertices, int numProcesses, int myRank) {
    assert(numProcesses >= 1 && myRank >= 0 && myRank < numProcesses);

    auto graph = allocateGraphPart(
            numVertices,
            getFirstGraphRowOfProcess(numVertices, numProcesses, myRank),
            getFirstGraphRowOfProcess(numVertices, numProcesses, myRank + 1)
    );

    if (graph == nullptr) {
        return nullptr;
    }

    assert(graph->numVertices > 0 && graph->numVertices == numVertices);
    assert(graph->firstRowIdxIncl >= 0 && graph->lastRowIdxExcl <= graph->numVertices);

    if (myRank == 0) {
        for (int hisRank = 0; hisRank < numProcesses; ++hisRank) {
            int firstRowIncl = getFirstGraphRowOfProcess(numVertices, numProcesses, hisRank);
            int lastRowExcl = getFirstGraphRowOfProcess(numVertices, numProcesses, hisRank + 1);
            for (int rowNum = firstRowIncl; rowNum < lastRowExcl; ++rowNum) {
                if (hisRank == 0) {
                    initializeGraphRow(graph->data[rowNum - firstRowIncl], rowNum, numVertices);
                } else {
                    initializeGraphRow(graph->extraRow, rowNum, numVertices);
                    MPI_Send(
                        graph->extraRow,
                        numVertices * sizeof(int),
                        MPI_BYTE,
                        hisRank,
                        0,
                        MPI_COMM_WORLD
                    );
                }
            }
        }
    } else {
        int firstRowIncl = getFirstGraphRowOfProcess(numVertices, numProcesses, myRank);
        int lastRowExcl = getFirstGraphRowOfProcess(numVertices, numProcesses, myRank + 1);
        for (int rowNum = firstRowIncl; rowNum < lastRowExcl; ++rowNum) {
            MPI_Recv(
                graph->data[rowNum - firstRowIncl],
                numVertices * sizeof(int),
                MPI_BYTE,
                0,
                0,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE
            );
        }
    }

    return graph;
}

void collectAndPrintGraph(Graph* graph, int numProcesses, int myRank) {
    assert(numProcesses >= 1 && myRank >= 0 && myRank < numProcesses);
    assert(graph->numVertices > 0);
    assert(graph->firstRowIdxIncl >= 0 && graph->lastRowIdxExcl <= graph->numVertices);

    if (myRank == 0) {
        for (int hisRank = 0; hisRank < numProcesses; ++hisRank) {
            int firstRowIncl = getFirstGraphRowOfProcess(graph->numVertices, numProcesses, hisRank);
            int lastRowExcl = getFirstGraphRowOfProcess(graph->numVertices, numProcesses, hisRank + 1);

            for (int numRow = firstRowIncl; numRow < lastRowExcl; ++numRow) {
                if (hisRank == 0) {
                    printGraphRow(graph->data[numRow], numRow, graph->numVertices);
                } else {
                    MPI_Recv(
                        graph->extraRow,
                        graph->numVertices * sizeof(int),
                        MPI_BYTE,
                        hisRank,
                        0,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                    );
                    printGraphRow(graph->extraRow, numRow, graph->numVertices);
                }
            }
        }
    } else {
        int firstRowIncl = getFirstGraphRowOfProcess(graph->numVertices, numProcesses, myRank);
        int lastRowExcl = getFirstGraphRowOfProcess(graph->numVertices, numProcesses, myRank + 1);
        for (int numRow = firstRowIncl; numRow < lastRowExcl; ++numRow) {
            MPI_Send(
                graph->data[numRow - firstRowIncl],
                graph->numVertices * sizeof(int),
                MPI_BYTE,
                0,
                0,
                MPI_COMM_WORLD
            );
        }
    }
}

void destroyGraph(Graph* graph, int numProcesses, int myRank) {
    freeGraphPart(graph);
}
