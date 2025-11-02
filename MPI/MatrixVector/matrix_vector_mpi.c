//
// Created by dominykas on 11/2/25.
//
#include "matrix_vector_mpi.h"

void dot_mpi(int argc, char **argv){
    //basic algo
    float* matrix_a;
    float* matrix_b;
    int n, m; //matrix size
    MPI_Init(&argc, &argv);
    int i = 0;
    int myRank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Barrier(MPI_COMM_WORLD);




}