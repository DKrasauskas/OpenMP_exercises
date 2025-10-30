//
// Created by domin on 20/10/2025.
//

#ifndef PARALLELIZED_SIMD_H
#define PARALLELIZED_SIMD_H

#include <omp.h>
#include <immintrin.h>
#include <fmaintrin.h>

#define SIMD_MULT(A, B) _mm256_mul_pd(A, B)
#define SIMD_FMA(A, B, C) _mm256_fmadd_pd(A, B, C)
#define SIMD_FDIV(A, B) _mm256_div_pd(A, B)
#define SIMD_FADD(A, B) _mm256_add_pd(A, B)
#define SIMD_STORE(A, B) _mm256_storeu_pd(A, B)


#define double4 __m256d
#define default_double4 _mm256_setzero_pd()
#define double4_SET(X) _mm256_set1_pd(X)
#define make_double4(A, B, C, D) _mm256_set_pd(A, B, C, D);

#endif //PARALLELIZED_SIMD_H
