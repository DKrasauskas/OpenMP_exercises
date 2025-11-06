//
// Created by dominykas on 11/5/25.
//

#ifndef OPENMP_EXERCISES_JACOBI_H
#define OPENMP_EXERCISES_JACOBI_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void swap_halos(int maxXCount, int maxYCount, double *src,
                int prevRank, int nextRank);


/*************************************************************
 * Performs one iteration of the Jacobi method and computes
 * the residual value.
 *
 * NOTE: u(0,*), u(maxXCount-1,*), u(*,0) and u(*,maxYCount-1)
 * are BOUNDARIES and therefore not part of the solution.
 *************************************************************/
double one_jacobi_iteration(double xStart, double yStart,
                            int maxXCount, int maxYCount,
                            double *src, double *dst,
                            double deltaX, double deltaY,
                            double alpha, double omega);
double one_jacobi_iteration_serial(double xStart, double yStart,
                            int maxXCount, int maxYCount,
                            double *src, double *dst,
                            double deltaX, double deltaY,
                            double alpha, double omega);
/**********************************************************
 * Checks the error between numerical and exact solutions
 **********************************************************/
double checkSolution(double xStart, double yStart, int maxXCount, int maxYCount, double *u, double deltaX, double deltaY, double apha);
double checkSolutionSerial(double xStart, double yStart,
                           int maxXCount, int maxYCount,
                           double *u,
                           double deltaX, double deltaY,
                           double alpha);
int jacobi(int argc, char **argv);
int jacobi_serial(int argc, char **argv);
#endif //OPENMP_EXERCISES_JACOBI_H
