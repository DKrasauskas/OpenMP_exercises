//
// Created by domin on 20/10/2025.
//
#include "pi.h"

__attribute__((always_inline)) double f(double a)
{
    return (4.0 / (1.0 + a*a));
}

int parallel_pi_generic(int thread_count, int n){
    double fSum[64][8];
    omp_set_num_threads(thread_count);
    int num_threads=thread_count;
    const double fH   = 1.0 / (double) n;
#pragma omp parallel shared(n, fSum) default(none)
    {
        const double fH   = 1.0 / (double) n;
        double TotalSum = 0.0;
        double fX;


        int tid = omp_get_thread_num();
        fSum[tid][0] = 0.0;
#pragma omp for schedule(static, n/20)
        for (int i = 0; i < n; i += 1)
        {
            fX = fH * ((double)i + 0.5);
            TotalSum += f(fX);
        }
        fSum[tid][0] = TotalSum;

    }
    double sum = 0;
    for (int j=0; j < num_threads; j += 1)
        sum = sum + fSum[j][0];
    printf("%f \n", sum * fH);
    return 0;
}

inline tuple generic_pi(int n){

    const double fH   = 1.0 / (double) n;
    double TotalSum = 0.0;
    double fSum;
    double fX;
    int i, j;

    fSum=0.0;
    for (i = 0; i < n; i += 1)
    {
        fX = fH * ((double)i + 0.5);
        fSum = fSum+ f(fX);
    }

    double pi = fH * TotalSum;
    tuple output;
    output.value = fSum * fH;
    return output;
}

inline tuple sequential_pi(int n){

        const double fH   = 1.0 / (double) n;
        int split = n;
        double sum = 0.0;
        int i;
        double4 vecSum = default_double4;
        double local_sum[4];

        for (i = 0; i < n - 5; i += 4)
        {
            double4 idx   = make_double4(i+3.5, i+2.5, i+1.5, i+0.5);
            double4 fX    = SIMD_MULT(idx,   double4_SET(fH));
            double4 denom = SIMD_FMA(fX, fX, double4_SET(1.0));
            double4 recip = SIMD_FDIV(double4_SET(4.0), denom);
            vecSum        = SIMD_FADD(vecSum, recip);

        }
        SIMD_STORE(local_sum, vecSum);
        sum += local_sum[0] + local_sum[1] + local_sum[2] + local_sum[3];
        double fX1;

        for (i; i <n; i += 1)
        {
            fX1 = fH * ((double)i + 0.5);
            double local_sum = (4.0 / (fma(fX1, fX1, 1.0)));
            sum += local_sum;
        }
    tuple output;
    output.value = sum * fH;
    return output;
}

inline tuple parallel_pi(int thread_count, int n){
    int num_threads = thread_count;
    omp_set_num_threads(num_threads);
    const double fH   = 1.0 / (double) n;
    double fSum[64][8];
    #pragma omp parallel shared(fSum, num_threads, n) default(none)
    {
        const double fH   = 1.0 / (double) n;
        int split = n / (int)num_threads;
        int tid = omp_get_thread_num();
        int begin = tid * split;
        int end = begin + split < n ? begin + split : n;
        double sum = 0.0;

        double4 vecSum = default_double4;
        double local_sum[4];
        int i = begin;

        for (i = begin; i <end - 4; i += 4)
        {
            double4 idx   = make_double4(i+3.5, i+2.5, i+1.5, i+0.5);
            double4 fX    = SIMD_MULT(idx,   double4_SET(fH));
            double4 denom = SIMD_FMA(fX, fX, double4_SET(1.0));
            double4 recip = SIMD_FDIV(double4_SET(4.0), denom);
            vecSum        = SIMD_FADD(vecSum, recip);

        }
        SIMD_STORE(local_sum, vecSum);
        sum += local_sum[0] + local_sum[1] + local_sum[2] + local_sum[3];
        double fX1;

        for (i; i <end; i += 1)
        {
            fX1 = fH * ((double)i + 0.5);
            double local_sum = (4.0 / (fma(fX1, fX1, 1.0)));
            sum += local_sum;
        }
        //only after the thread is done shall it write to the shared memory:
        fSum[tid][0] = sum;

    }
    double sum = 0.0;
    for(int i = 0; i < num_threads; i++){
        sum += fSum[i][0];
    }
    tuple output;
    output.value = sum * fH;
    return output;
}

//    int num_threads = 20;
//    omp_set_num_threads(num_threads);
//    printf("max %d\n", omp_get_max_threads());
//    double fSum[64][8];
//
//    int n = 1000000000;
//    const double fH0   = 1.0 / (double) n;
//    double fTimeStart, fTimeEnd;
//    fTimeStart = omp_get_wtime();
//
//    #pragma omp parallel shared(fSum, num_threads, n) default(none)
//    {
//        const double fH   = 1.0 / (double) n;
//        int split = n / (int)num_threads;
//        int tid = omp_get_thread_num();
//        int begin = tid * split;
//        int end = begin + split < n - 1 ? begin + split : n - 1;
//        double sum = 0.0;
//
//        __m256d vecSum = _mm256_setzero_pd();
//        double local_sum[4];
//        int i = begin;
//        #pragma unroll
//        for (i = begin; i <end - 4; i += 4)
//        {
//            __m256d idx = _mm256_set_pd(i+3.5, i+2.5, i+1.5, i+0.5);
//            __m256d fX = _mm256_mul_pd(idx, _mm256_set1_pd(fH));
//            __m256d denom = _mm256_fmadd_pd(fX, fX, _mm256_set1_pd(1.0));
//            __m256d recip = _mm256_div_pd(_mm256_set1_pd(4.0), denom);
//            vecSum = _mm256_add_pd(vecSum, recip);
//
//        }
//        _mm256_storeu_pd(local_sum, vecSum);
//        sum += local_sum[0] + local_sum[1] + local_sum[2] + local_sum[3];
//        double fX1;
//
//        for (i; i <end; i += 1)
//        {
//            fX1 = fH * ((double)i + 0.5);
//            double local_sum = (4.0 / (fma(fX1, fX1, 1.0)));
//            sum += local_sum;
//        }
//        //only after the thread is done shall it write to the shared memory:
//        fSum[tid][0] = sum;
//
//    }
//    double sum = 0.0;
//    for(int i = 0; i < num_threads; i++){
//        sum += fSum[i][0];
//    }
//    fTimeEnd = omp_get_wtime();
//    printf("  wall clock time     = %.20f\n", fTimeEnd - fTimeStart);
//    printf("%.20f", sum * fH0);
//    printf("Hello, World!\n");
#include "pi.h"
