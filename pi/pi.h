//
// Created by domin on 20/10/2025.
//

#ifndef PARALLELIZED_PI_H
#define PARALLELIZED_PI_H
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include "SIMD.h"
#include "../typedefs.h"


//const double fPi25DT = 3.141592653589793238462643;

__attribute__((always_inline)) double f(double a);

int parallel_pi_generic(int thread_count, int n);
tuple sequential_pi(int n);
tuple generic_pi(int n);
tuple parallel_pi(int thread_count, int n);

#endif //PARALLELIZED_PI_H
