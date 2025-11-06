//
// Created by dominykas on 11/5/25.
//

#ifndef OPENMP_EXERCISES_MATRIX_MATRIX_H
#define OPENMP_EXERCISES_MATRIX_MATRIX_H
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
void mat_mat();
void mat_mat_benchmark(float* a, float*b, float* c, int ns, int thread_num);
void mat_mat_parallel();
#endif //OPENMP_EXERCISES_MATRIX_MATRIX_H
