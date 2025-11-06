//
// Created by domin on 30/10/2025.
//
#include "benchmarking.h"
#include "openMP/pi/pi.h"

void benchmark(tuple (*func)(int),  int n){
    double begin = omp_get_wtime();
    tuple output = (*func)(n);
    double end = omp_get_wtime();
    double dt = end - begin;
    printf("%.20f\n", dt);
}

void benchmarki(tuple (*func)(int, int), int threads, int n){
    double begin = omp_get_wtime();
    tuple output =(*func)(threads, n);
    double end = omp_get_wtime();
    double dt = end - begin;
    printf("%.20f\n", dt);
}

void benchmarkv(tuple (*func)(int, int), int threads, int n){
    double begin = omp_get_wtime();
    tuple output =(*func)(threads, n);
    double end = omp_get_wtime();
    double dt = end - begin;
    printf("%.20ff\n", dt);
}