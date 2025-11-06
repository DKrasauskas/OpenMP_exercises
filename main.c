
#include <stdio.h>
#include <omp.h>
#include "openMP/pi/pi.h"
#include "benchmarking.h"
#include "openMP/fixit/fixit.h"
#include "openMP/matrix_vector/mv.h"
#include "MPI/integration/integration.h"
#include "utils/utils.h"
#include "MPI/ping_pong/pingPong-a.h"
#include "MPI/Jacobi/jacobi.h"
#include "openMP/matrix_matrix/matrix-matrix.h"
#include "MPI/MatrixVector/matrix_vector_mpi.h"


void get_runtime_mpi_matrix_matrix(int argc, char** argv){
    MPI_Init(&argc, &argv);
    int myRank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Barrier(MPI_COMM_WORLD);
    matrix3 mat;
    mat.n = 1000;
    mat.m = 1000;
    mat.p = 1000;
    if(myRank == 0){
        mat.matrix_a = (float*)malloc(sizeof(float) * mat.m * mat. n);
        mat.matrix_b = (float*)malloc(sizeof(float) * mat.m * mat. p);
        mat.matrix_c = (float*)malloc(sizeof(float) * mat.p * mat. n);

    }
    MPI_Barrier(MPI_COMM_WORLD);
    double begin = MPI_Wtime();
    dot_mpi(mat);
    double end = MPI_Wtime();
    if(myRank == 0){
        printf("took %f s \n", end - begin);
        free(mat.matrix_a);
        free(mat.matrix_b);
        free(mat.matrix_c);
    }

    MPI_Finalize();
}

void get_runtime_openmp_matrix_matrix(int argc, char** argv) {
    matrix3 mat;
    mat.n = 1000;
    mat.m = 1000;
    mat.p = 1000;
    mat.matrix_a = (float*)malloc(sizeof(float) * mat.m * mat. n);
    mat.matrix_b = (float*)malloc(sizeof(float) * mat.m * mat. p);
    mat.matrix_c = (float*)malloc(sizeof(float) * mat.p * mat. n);
    for(int i = 1; i < 20; i += 2){
        double begin = omp_get_wtime();
        mat_mat_benchmark(mat.matrix_a, mat.matrix_b, mat.matrix_c, mat.n,  20);
        double end = omp_get_wtime();
        printf("%.10f \n", end - begin);
    }

}
int main(int argc, char** argv) {
//    int p[256];
//    printf("%d", p);
//    int n = 1500000000;
//    int max_threads = omp_get_max_threads();
//    printf("runtime generic sequential:\n");
//    benchmark(generic_pi,  n);
//    #ifdef MVX_ARCH
//    benchmark(sequential_pi,  n);
//    printf("runtime parallel:\n");
//    #endif
//    for(int i = 1; i < max_threads + 1; i ++){
//        benchmarki(reduction_parallel_pi, i, n);
//    }
//    return 0;
    generate_ab(argc, argv);
}
//13.34123113
//5.268928s
//1.653267s
//1.292311