#include <cmath>
#include<immintrin.h>
#include<iostream>

namespace reference {
    /*
    Flop count:
    additions: 2 * blockDimension²
    multiplication: 1 * blockDimension² + 2
    sqrt: 1
    Total: 3 * blockDimension² + 3
    */

    // used by base implementation (compress.cpp)
    void difference_norm(const int blockDimension, double *ar1, double *ar2, double *sum) {
        *sum = 0;
        for(int i = 0; i < blockDimension * blockDimension; i++){
            *sum += std::pow(ar1[i] - ar2[i], 2);
        }
        *sum = std::sqrt((*sum) / (blockDimension*blockDimension)) / 255.0;
    }
    
    void difference_norm(const int blockDimension, float *ar1, float *ar2, float *sum) {
        *sum = 0.f;
        for(int i = 0; i < blockDimension * blockDimension; i++){
            *sum += powf(ar1[i] - ar2[i], 2);
        }
        *sum = sqrtf((*sum) / (blockDimension*blockDimension)) / 255.0f;
    }
}


namespace fast {
    void difference_norm(const int blockDimension, double *ar1, double *ar2, double *sum) {
        *sum = 0.;
        for(int i = 0; i < blockDimension * blockDimension; i++){
            *sum += std::pow(ar1[i] - ar2[i], 2);
        }
        *sum = std::sqrt((*sum) / (blockDimension*blockDimension)) / 255.0;
    }
    
    void difference_norm(const int blockDimension, float *ar1, float *ar2, float *sum) {
        *sum = 0.f;
        for(int i = 0; i < blockDimension * blockDimension; i++){
            *sum += powf(ar1[i] - ar2[i], 2);
        }
        *sum = sqrtf((*sum) / (blockDimension*blockDimension)) / 255.0f;
    }

    void difference_norm_scalar_replacement(const int blockDimension, float *ar1, float *ar2, float *sum) {
        float s = 0.f;
        for(int i = 0; i < blockDimension * blockDimension; i++){
            s += powf(ar1[i] - ar2[i], 2);
        }
        *sum = sqrtf((s) / (blockDimension*blockDimension)) / 255.0f;
    }
    
    void difference_norm_scalar_replacement_modified(const int blockDimension, float *ar1, float *ar2, float *sum) {
        //idea: does not actually compute the same error.
        //but we are only interested inf finding the smallest error and therefore can omit:
        //division by numbers >= 1 (monotone)
        //sqrt (monotone)
        float s = 0.f;
        for(int i = 0; i < blockDimension * blockDimension; i++){
            const float x1 = ar1[i];
            const float x2 = ar2[i];
            const float difference = x1 - x2;
            const float difference_squared = difference * difference;
            s += difference_squared;
        }
        *sum = s;
    }

    void difference_norm_vectorized(const int blockDimension, float *ar1, float *ar2, float *sum) {
        float s = 0.f;
        __m256 s_vec,
        ar1_vec, ar2_vec,
        temp;

        s_vec = _mm256_setzero_ps();

        for(int i = 0; i < blockDimension * blockDimension; i+=8){
            ar1_vec = _mm256_load_ps(ar1 + i);
            ar2_vec = _mm256_load_ps(ar2 + i);
            temp = _mm256_sub_ps(ar1_vec, ar2_vec);
            temp = _mm256_mul_ps(temp, temp);
            s_vec = _mm256_add_ps(s_vec, temp);
        }
        s_vec = _mm256_hadd_ps(s_vec, _mm256_setzero_ps());
        s_vec = _mm256_hadd_ps(s_vec, _mm256_setzero_ps());
        temp = _mm256_permute2f128_ps(s_vec, _mm256_setzero_ps(),_MM_SHUFFLE(0,0,0,1));
        s_vec = _mm256_add_ps(s_vec, temp);

        s = _mm256_cvtss_f32(s_vec);
        *sum = sqrtf((s) / (blockDimension*blockDimension)) / 255.0f;
    }
    
    void difference_norm_modified(const int blockDimension, float *ar1, float *ar2, float *sum) {
        //idea: does not actually compute the same error.
        //but we are only interested inf finding the smallest error and therefore can omit:
        //division by numbers >= 1 (monotone)
        //sqrt (monotone)
        float s = 0.f;
        #pragma GCC ivdep
        for(int i = 0; i < blockDimension * blockDimension; i++){
            const float x1 = ar1[i];
            const float x2 = ar2[i];
            const float difference = x1 - x2;
            const float difference_squared = difference * difference;
            s += difference_squared;
        }
        *sum = s;
    }

    void difference_norm_vectorized_without_normalization(const int blockDimension, float *ar1, float *ar2, float *sum) {
        float s = 0.f;
        __m256 s_vec,
        ar1_vec, ar2_vec,
        temp;

        s_vec = _mm256_setzero_ps();

        for(int i = 0; i < blockDimension * blockDimension; i+=8){
            ar1_vec = _mm256_load_ps(ar1 + i);
            ar2_vec = _mm256_load_ps(ar2 + i);
            temp = _mm256_sub_ps(ar1_vec, ar2_vec);
            temp = _mm256_mul_ps(temp, temp);
            s_vec = _mm256_add_ps(s_vec, temp);
        }
        s_vec = _mm256_hadd_ps(s_vec, _mm256_setzero_ps());
        s_vec = _mm256_hadd_ps(s_vec, _mm256_setzero_ps());
        temp = _mm256_permute2f128_ps(s_vec, _mm256_setzero_ps(),_MM_SHUFFLE(0,0,0,1));
        s_vec = _mm256_add_ps(s_vec, temp);

        s = _mm256_cvtss_f32(s_vec);
        *sum = s;
    }
}
