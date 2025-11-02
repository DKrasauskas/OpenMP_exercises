//
// Created by dominykas on 11/2/25.
//

#ifndef OPENMP_EXERCISES_INTEGRATION_H
#define OPENMP_EXERCISES_INTEGRATION_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define IDLETIME 0.1

#define TAG_WORK 0
#define TAG_END 2

double func(double x);

double integrate(double (*f)(double x), double x_start, double x_end, int maxSteps);

int integration_dispatch(int argc, char **argv);
#endif //OPENMP_EXERCISES_INTEGRATION_H
