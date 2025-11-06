//
// Created by dominykas on 11/2/25.
//
#include "matrix_vector_mpi.h"


matrix3 dot_mpi(matrix3 matrices){
    //basic algo
    float* matrix_a = matrices.matrix_a;
    float* matrix_c =  matrices.matrix_c;
    int n =  matrices.n;
    int m =  matrices.m;
    int p =  matrices.p;
    int myRank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Barrier(MPI_COMM_WORLD);
    int local_rows = n / numProcs;
    int global_row_index;
    float* subarray = (float*)malloc(sizeof(float) * local_rows * m);
    float* result   = (float*)malloc(sizeof(float) * local_rows * p);
    float* matrix_b =  myRank == 0 ? matrices.matrix_b :  (float*)malloc(sizeof(float) * p * m);
    MPI_Barrier(MPI_COMM_WORLD);

    //now, time to broadcast the matrix b to the other processes:
    MPI_Bcast(matrix_b, p * m, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //now scatter the initial array
    MPI_Scatter(matrix_a, local_rows * m, MPI_FLOAT, subarray, local_rows * m, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //now its time to perform matrix - matrix operations for the sub-processes

    for(int i = 0; i < local_rows; i ++){
        //iterate through columns of matrix b
        for(int j = 0; j < p; j ++){
            float dotp = 0.0f;
            for(int k = 0; k < m; k ++){
                int index_b = j + k * p;
                int index_a = k + i * m;
                dotp += matrix_b[index_b] * subarray[index_a];
            }
            //load to the resulting
            int index_result = j + i * p;
            result[index_result] = dotp;
        }
    }
    //now gather all the results:
    MPI_Gather(result, p * local_rows, MPI_FLOAT, matrix_c, p * local_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);
    free(subarray);
    free(result);
    if(myRank != 0){
        free(matrix_b);
    }
    return matrices;
}


matrix3 dot_mpi_sequential(matrix3 matrices){
    float* matrix_a  = matrices.matrix_a;
    float* matrix_b  = matrices.matrix_b;
    float* matrix_c =  matrices.matrix_c;
    int n =  matrices.n;
    int m =  matrices.m;
    int p =  matrices.p;
    for(int i = 0; i < n; i ++){
        //iterate through columns of matrix b
        for(int j = 0; j < p; j ++){
            float dotp = 0.0f;
            for(int k = 0; k < m; k ++){
                int index_b = j + k * p;
                int index_a = k + i * m;
                dotp += matrix_b[index_b] * matrix_a[index_a];
            }
            //load to the resulting
            int index_result = j + i * p;
            matrix_c[index_result] = dotp;
        }
    }
    return matrices;
}