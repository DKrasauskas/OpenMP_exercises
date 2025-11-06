//
// Created by dominykas on 11/2/25.
//
/******************************************************************************
* FILE: fixit.c
*
*   This very simple program contains several errors. Find them and fix.
*   For instance, is there dependency?
*
******************************************************************************/
#include "fixit.h"


#define __syncthreads() _Pragma("omp barrier")

void fixit()
{
    int nrthreads, tid, i, j;
    //might run out of stack
    int a[N][N], b[N][N];
    //have the master thread do the initialization BEFORE private assign
    for (j = 0; j < N; j++) a[0][j] = j;
    /* Fork a team of threads */
    #define __syncthreads() _Pragma("omp barrier")
    #pragma omp parallel shared(nrthreads) private(i, a, tid, j, b)
    {
        /* Obtain/print thread info */
        /// tid needs to be private
        tid = omp_get_thread_num();
        if (tid == 0) {
            nrthreads = omp_get_num_threads();
            printf("Number of threads = %d\n", nrthreads);
            for (j = 0; j < N; j++) a[0][j] = j;
        }
        __syncthreads();
        for (j = 0; j < N; j++) a[0][j] = j;

        for (i = 1; i < N; i++) {
            for (j = 0; j < N; j++) {
                a[i][j] = nrthreads * i + j;
                b[i][j] = a[i - 1][j] + a[i][j];
            }
        }



       printf("Thread %d done. Last element of a[%d][%d]= %d\n",tid,tid*N/nrthreads, N-1, a[tid*N/nrthreads][N-1]);
       // printf("Thread %d done. Values of b[%d][%d]= %d\n",tid,tid*N/nrthreads, N-1, b[tid*N/nrthreads][N-1]);


     /* All threads join master thread and disband */
    }
}



