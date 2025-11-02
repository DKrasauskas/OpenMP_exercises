//
// Created by dominykas on 11/2/25.
//

#ifndef OPENMP_EXERCISES_UTILS_H
#define OPENMP_EXERCISES_UTILS_H
/*
 *  OpenMP Exercise
 *
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>

double CalcPi(int n, int threads);

int call()
{
    int n = 15E7;  /* n is large number 150000000  */
    const double fPi25DT = 3.141592653589793238462643;
    double fPi;
    double fTimeStart, fTimeEnd;

#ifdef READ_INPUT
    printf("Enter the number of intervals: ");
    scanf("%d",&n);
#endif

    if (n <= 0 || n > 2147483647 )
    {
        printf("\ngiven value has to be between 0 and 2147483647\n");
        return 1;
    }

    int threads_max = omp_get_max_threads();
    for (int i = 1; i < threads_max; i++)
    {
        fTimeStart = omp_get_wtime();

        fPi = CalcPi(n, i);

        fTimeEnd = omp_get_wtime();
        printf("%.20f\n", fTimeEnd - fTimeStart);
    }
    return 0;
}


__attribute__((always_inline)) double fa(double a)
{
    return (4.0 / (1.0 + a*a));
}


double CalcPi(int n, int threads)
{
    const double fH   = 1.0 / (double) n;
    int i, j, num_threads ;
    double sum = 0.0;

    omp_set_num_threads(threads);

//start of parallel region

    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < n; i += 1)
    {
        double fX = fH * ((double)i + 0.5);
        sum += fa(fX);
    }

//end of parallel region
    return (fH * sum);
}







#endif //OPENMP_EXERCISES_UTILS_H
