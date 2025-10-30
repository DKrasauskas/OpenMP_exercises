//
// Created by domin on 30/10/2025.
//
#include "benchmarking.h"
void benchmark(tuple (*func)(int),  int n){
    double begin = omp_get_wtime();
    tuple output = (*func)(n);
    double end = omp_get_wtime();
    double dt = end - begin;
    //printf("%.20f'\n", (output.value - fPi25DT) / fPi25DT * 100);
    printf("%.20f\n", dt);
}

void benchmarki(tuple (*func)(int, int), int threads, int n){
    double begin = omp_get_wtime();
    tuple output =(*func)(threads, n);
    double end = omp_get_wtime();
    double dt = end - begin;
    // printf("%.20f'\n", (output.value - fPi25DT) / fPi25DT * 100);
    printf("%.20f\n", dt);
}

void benchmarkv(tuple (*func)(int, int), int threads, int n){
    double begin = omp_get_wtime();
    tuple output =(*func)(threads, n);
    double end = omp_get_wtime();
    double dt = end - begin;
    printf("%.20ff\n", dt);
}