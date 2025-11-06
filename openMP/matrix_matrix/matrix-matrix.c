//
// Created by dominykas on 11/5/25.
//
/******************************************************************************
* FILE: mm.c
* DESCRIPTION:
*   This program calculates the product of matrix a[nra][nca] and b[nca][ncb],
*   the result is stored in matrix c[nra][ncb].
*   The max dimension of the matrix is constraint with static array declaration,
*   for a larger matrix you may consider dynamic allocation of the arrays,
*   but it makes a parallel code much more complicated (think of communication),
*   so this is only optional.
*
******************************************************************************/

#include "matrix-matrix.h"

// #define NRA 896                 /* number of rows in matrix A */
// #define NCA 896                 /* number of columns in matrix A */
// #define NCB 896                  /* number of columns in matrix B */


#define FOR _Pragma("omp parallel for")

void mat_mat_benchmark(float* a, float*b, float* c, int ns, int thread_num)
{
    float* matrix_a  = a;
    float* matrix_b  = b;
    float* matrix_c =  c;
    int n =  ns;
    int m =  ns;
    int p =  ns;
    omp_set_num_threads(thread_num);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < n; i ++){
        //iterate through columns of matrix b
        for(int j = 0; j < p; j ++){
            float dotp = 0.0f;
            for(int k = 0; k < m; k ++){
                int index_b = j + k * p;
                int index_a = k + i * m;
                dotp += matrix_b[index_b] * matrix_a[index_a];
            }
            //load to the resulting
            int index_result = j + i * p;
            matrix_c[index_result] = dotp;
        }
    }
}

void mat_mat()
{
    for(int N = 0; N < 400; N += 10) {
        int tid, nthreads, i, j, k, n, chunk;
/* for simplicity, set NRA=NCA=NCB=N  */
        int NRA = N;
        int NCA = N;
        int NCB = N;

        double a[NRA][NCA],           /* matrix A to be multiplied */
        b[NCA][NCB],           /* matrix B to be multiplied */
        c[NRA][NCB];           /* result matrix C */

        double fTimeStart, fTimeEnd, CalcTime;

/* set loop iteration chunk size, experiment with several different values of chunk */
        chunk = 8;

        /*** Initialize matrices ***/

        for (i = 0; i < NRA; i++)
            for (j = 0; j < NCA; j++)
                a[i][j] = i + j;

        for (i = 0; i < NCA; i++)
            for (j = 0; j < NCB; j++)
                b[i][j] = i * j;

        for (i = 0; i < NRA; i++)
            for (j = 0; j < NCB; j++)
                c[i][j] = 0;

        /* Use OpenMP to parallelize the computation of the following matrix-matrix multiplication. */
        //get sequential results
        fTimeStart = omp_get_wtime();

        /* for small N, the execution time can vary quite a lot, so we repeat the matrix-product computation 50 times */
        if (N < 4800)
            for (n = 0; n < 50; n++) {
                for (i = 0; i < NRA; i++) {
                    for (j = 0; j < NCB; j++)
                        for (k = 0; k < NCA; k++)
                            c[i][j] += a[i][k] * b[k][j];
                }
            }  /* end of if for-loop n */
        else  /* N>=48, we do matrix-matrix multiplication only once */
        {
            /* To Do: Pararallelize the computation of the matrix-matrix multiplication. */
            for (i = 0; i < NRA; i++) {
                for (j = 0; j < NCB; j++)
                    for (k = 0; k < NCA; k++)
                        c[i][j] += a[i][k] * b[k][j];
            }
        }

        fTimeEnd = omp_get_wtime();
        if (N < 48) {
            CalcTime = (fTimeEnd - fTimeStart) / 50;
        } else {
            CalcTime = fTimeEnd - fTimeStart;
        }
        printf(" %.20f\n", CalcTime);
//        printf("  c[0][0]    = %.20f\n", c[0][0]);
//        printf("  c[0][10]    = %.20f\n", c[0][10]);
//        printf("  c[100][0]    = %.20f\n", N < 100 ? c[N - 1][10] : c[100][10]);
    }
}


void mat_mat_parallel()
{
    for(int N = 0; N < 400; N += 10){
        int tid, nthreads, i, j, k, n, chunk;
/* for simplicity, set NRA=NCA=NCB=N  */
        int NRA = N;
        int NCA = N;
        int NCB = N;

        double a[NRA][NCA],           /* matrix A to be multiplied */
        b[NCA][NCB],           /* matrix B to be multiplied */
        c[NRA][NCB];           /* result matrix C */

        double fTimeStart, fTimeEnd, CalcTime;

/* set loop iteration chunk size, experiment with several different values of chunk */
        chunk = 8;

        /*** Initialize matrices ***/

        for (i = 0; i < NRA; i++)
            for (j = 0; j < NCA; j++)
                a[i][j] = i + j;

        for (i = 0; i < NCA; i++)
            for (j = 0; j < NCB; j++)
                b[i][j] = i * j;

        for (i = 0; i < NRA; i++)
            for (j = 0; j < NCB; j++)
                c[i][j] = 0;

        /* Use OpenMP to parallelize the computation of the following matrix-matrix multiplication. */
        //get sequential results
        omp_set_num_threads(4);
        fTimeStart = omp_get_wtime();

        /* for small N, the execution time can vary quite a lot, so we repeat the matrix-product computation 50 times */
        if (N < 1000)
            for (n = 0; n < 50; n++) {
                #pragma omp parallel for schedule(dynamic)
                for (int i = 0; i < NRA; i++) {
                    //subdivide rows among processors and then compute the dot products with matrix B
                    for (int k = 0; k < NCA; k++) {
                        for (int j = 0; j < NCB; j++) {
                            c[i][j] += a[i][k] * b[k][j];
                        }
                    }
                }
            }  /* end of if for-loop n */
        else  /* N>=48, we do matrix-matrix multiplication only once */
        {
        #pragma omp parallel for schedule(static, chunk)
            for (int i = 0; i < NRA; i++) {
                //subdivide rows among processors and then compute the dot products with matrix B
                for (int k = 0; k < NCA; k++) {
                    for (int j = 0; j < NCB; j++) {
                        c[i][j] += a[i][k] * b[k][j];
                    }
                }
            }
        }

        fTimeEnd = omp_get_wtime();
        if (N < 48) {
            CalcTime = (fTimeEnd - fTimeStart) / 50;
        } else {
            CalcTime = fTimeEnd - fTimeStart;
        }
        printf("%.20f\n", CalcTime);
//        printf("  c[0][0]    = %.20f\n", c[0][0]);
//        printf("  c[0][10]    = %.20f\n", c[0][10]);
        //printf("  c[100][0]    = %.20f\n", N < 100 ? c[N - 1][10] : c[100][10]);
    }
}