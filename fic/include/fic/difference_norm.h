#if !defined(__FIC_DIFFERENCE_NORM_H)
#define __FIC_DIFFERENCE_NORM_H

namespace reference {
    void difference_norm(const int, double *, double *, double *);
    void difference_norm(const int, float *, float *, float *);
}

namespace fast {
    void difference_norm(const int, double *, double *, double *);
    void difference_norm(const int, float *, float *, float *);
    void difference_norm_scalar_replacement(const int, float *, float *, float *);
    void difference_norm_scalar_replacement_modified(const int, float *, float *, float *);
    void difference_norm_vectorized(const int, float *, float *, float *);
    void difference_norm_vectorized_without_normalization(const int, float *, float *, float *);
    void difference_norm_vectorized_modified(const int, float *, float *, float *);
    void difference_norm_vectorized_without_normalization(const int, float *, float *, float *);
}

#endif // __FIC_DIFFERENCE_NORM_H
