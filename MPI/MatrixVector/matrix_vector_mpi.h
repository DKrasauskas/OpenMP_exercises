//
// Created by dominykas on 11/2/25.
//

#ifndef OPENMP_EXERCISES_MATRIX_VECTOR_MPI_H
#define OPENMP_EXERCISES_MATRIX_VECTOR_MPI_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

typedef struct{
    float* matrix_a, *matrix_b, *matrix_c;
    int n,  m,  p;
}matrix3;

matrix3 dot_mpi(matrix3 matrices);
matrix3 dot_mpi_sequential(matrix3 matrices);
#endif //OPENMP_EXERCISES_MATRIX_VECTOR_MPI_H
