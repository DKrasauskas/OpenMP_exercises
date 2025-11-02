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
#include <string.h>
#include "mv.h"
#include "../pi/SIMD.h"

void init(int i, int j, int size, float* matrix, float* b, float* c){
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            int index = j + i * size;
            matrix[index] = fminf(i * 1.0 / (j + 1.0), j * 1.0 / (i + 1.0));
        }
        b[i] = 1.0 * (i + 1);
        c[i] = 0.0;
    }
}
void wrong_parallelization(int i, int j, int size, float* matrix, float* b, float* c){
    #pragma omp parallel for private(i)
    for (i = 0; i < SIZE; i++){
        for (j = 0; j < SIZE; j++) {
            int index = i + j * size;
            c[j] = c[j] + matrix[index] * b[i];
        }
    }
}

void correct_parallelization(int size, float* matrix, float* b, float* c){
        #pragma omp parallel for
        for (int i=0; i < SIZE; i++){
            for (int j=0; j < SIZE; j++){
                int index = j + i * size;
                c[i] = c[i] + matrix[index] * b[j];
            }
        }
}

void fast_parallelization(int size, float* matrix, float* b, float* c){
    int threads = omp_get_max_threads();
    int split = size / threads;
    float8 *restrict matrix_cast  = (float8*)matrix;
    float8 *restrict b_cast =  (float8*)b;
    #pragma omp parallel private(threads)
    {
        int tid = omp_get_thread_num();
        //each processor gets n / p rows:
        int begin = split * tid; //begin row
        int end = tid == omp_get_max_threads() - 1 ? size : begin + split; //end row

        for(int i = begin; i < end; i ++){
            float8 result = default_float8;
            int j;
            for(j = 0; j < size - 8; j += 8){
                int index = (j + i * size) / 8;
                float8 elements = matrix_cast[index];
                float8 b_vectorized = b_cast[j / 8];
                result = SIMD_FMA_PS(elements, b_vectorized, result);
            }
            //now gather the results:
            float temp[8];
            SIMD_STORE_PS(temp, result);
            float sum = 0;
            #pragma unroll
            for(int p = 0; p < 8; p ++) sum += temp[p];
            //now compute from the end
            for(int k = j; k < size; k ++){
                int index = j + i * size;
                sum += matrix[index] * b[k];
            }
            c[i] = sum;
        }
    }
}

int dot () {
    int size = SIZE;
    int i, j;
    double fTimeStart, fTimeEnd;

    //i have limited stack so dynamic alloc is the way to go
    int padded_size = ((size * size)/ 8) * 8 + 8;
    int padded_b_size = ((size)/ 8) * 8 + 8;
    float *matrix;
    posix_memalign((void**)&matrix, 32, padded_size * sizeof(float));
    float *b = (float *) malloc(sizeof(float) * size);
    posix_memalign((void**)&b, 32, padded_b_size * sizeof(float));
    float *c = (float *) malloc(sizeof(float) * size );
    float *target = (float *) malloc(sizeof(float) * size );

    /* Initializations */
    init(i, j, size,  matrix, b, c);

    fTimeStart = omp_get_wtime();

    //#pragma omp parallel for private(i)
    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++) {
            int index = j + i * size;
            c[j] = c[j] + matrix[index] * b[i];
        }
    }

    fTimeEnd = omp_get_wtime();
    memcpy(target, c, sizeof(float) * size);
    printf("  wall clock time     = %.20f\n", fTimeEnd - fTimeStart);
    printf("runtime slow algorithm\n");
    int k = omp_get_max_threads();
    for(int t = 1; t < k; t ++){
        omp_set_num_threads(t);
        init(i, j, size,  matrix, b, c);
        fTimeStart = omp_get_wtime();
        correct_parallelization(size, matrix, b, c);
        fTimeEnd = omp_get_wtime();
        printf("%.20f\n", fTimeEnd - fTimeStart);
        //find absolute error
//        float max = -0xFFFF;
//        for (j = 0; j < size; j++) {
//           float local_err = (c[j] - target[j]) / target[j] * 100;
//           if(local_err > max) max = local_err;
//        }
//        printf("max err = %f \n", max);
    }

    printf("runtime fast algorithm:\n");
    fTimeStart = omp_get_wtime();
    for(int t = 1; t < k; t ++){
        omp_set_num_threads(t);
        init(i, j, size,  matrix, b, c);
        fTimeStart = omp_get_wtime();
        fast_parallelization(size, matrix, b, c);
        fTimeEnd = omp_get_wtime();
        printf("%.20f\n", fTimeEnd - fTimeStart);
    }

    fTimeEnd = omp_get_wtime();

    printf("\n");
    printf("  wall clock time     = %.20f\n", fTimeEnd - fTimeStart);
    free(matrix);
    free(c);
    free(b);
    free(target);
}


