//
// Created by domin on 30/10/2025.
//

#ifndef OPENMP_EXERCISES_BENCHMARKING_H
#define OPENMP_EXERCISES_BENCHMARKING_H
#include "typedefs.h"
#include <omp.h>
#include <stdio.h>

void benchmark(tuple (*func)(int),  int n);

void benchmarki(tuple (*func)(int, int), int threads, int n);

void benchmarkv(tuple (*func)(int, int), int threads, int n);

#endif //OPENMP_EXERCISES_BENCHMARKING_H
