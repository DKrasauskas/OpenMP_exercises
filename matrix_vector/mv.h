//
// Created by dominykas on 11/2/25.
//

#ifndef OPENMP_EXERCISES_MV_H
#define OPENMP_EXERCISES_MV_H
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#define SIZE 1000

//
// Created by dominykas on 11/2/25.
//
/******************************************************************************
* Parallel Matrix-Vector multiplication - C/C++ Version
* FILE: Parallel-MV.c
* DESCRIPTION:
*   This example tries to parallelize the multiplication of matrix A with a vector, and
*   stores the result in vector c.
******************************************************************************/
#include "mv.h"
#include "../pi/SIMD.h"
void wrong_parallelization(int i, int j, int size, float* matrix, float* b, float* c);

void correct_parallelization(int size, float* matrix, float* b, float* c);

void fast_parallelization(int size, float* matrix, float* b, float* c);

int dot ();




#endif //OPENMP_EXERCISES_MV_H
