//
// Created by domin on 20/10/2025.
//

#ifndef PARALLELIZED_SIMD_H
#define PARALLELIZED_SIMD_H

#include <omp.h>
#include <immintrin.h>
#include <fmaintrin.h>
#define MVX_ARCH
#ifdef MVX_ARCH
#define SIMD_MULT(A, B) _mm256_mul_pd(A, B)
#define SIMD_FMA(A, B, C) _mm256_fmadd_pd(A, B, C)
#define SIMD_FDIV(A, B) _mm256_div_pd(A, B)
#define SIMD_FADD(A, B) _mm256_add_pd(A, B)
#define SIMD_STORE(A, B) _mm256_storeu_pd(A, B)


#define double4 __m256d
#define default_double4 _mm256_setzero_pd()
#define double4_SET(X) _mm256_set1_pd(X)
#define make_double4(A, B, C, D) _mm256_set_pd(A, B, C, D);
//and also for floats
#define SIMD_MULT_PS(A, B) _mm256_mul_ps(A, B)
#define SIMD_FMA_PS(A, B, C) _mm256_fmadd_ps(A, B, C)
#define SIMD_FDIV_PS(A, B) _mm256_div_ps(A, B)
#define SIMD_FADD_PS(A, B) _mm256_add_ps(A, B)
#define SIMD_STORE_PS(A, B) _mm256_storeu_ps(A, B)

#define float8 __m256
#define default_float8 _mm256_setzero_ps()
#define float8_SET(X) _mm256_set1_ps(X)
#define make_float8(A, B, C, D, E, F, G, H) _mm256_set_ps(A, B, C, D, E, F, G, H)
#endif
#endif //PARALLELIZED_SIMD_H
