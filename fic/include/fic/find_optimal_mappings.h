#if !defined(__FIC_FIND_OPTIMAL_MAPPINGS_H)
#define __FIC_FIND_OPTIMAL_MAPPINGS_H

namespace reference {
    void find_optimal_mappings(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_flip(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_with_entropy(const int, const int, const int, const int, float *, float **, float *);
}

namespace fast {
    void find_optimal_mappings(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_scalar_replacement(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_scalar_replacement_modified(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_inline_modified(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_strength_reduction_modified(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_code_motion_modified(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_ilp_modified(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_ilp_modified_reordered(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_GCC_vectorization_modified(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_avx(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_avx_splitloop(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_avx_althsum(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_avx_hsuminlined(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_avx_reordered(const int, const int, const int, const int, float *, float **, float *);
    

    void find_optimal_mappings_flip(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_flip_scalar_replacement_modified(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_flip_ilp_modified(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_flip_ilp_modified_splitloop(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_flip_ilp_modified_splitloop2(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_flip_ilp_modified_reordered(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_flip_avx(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_flip_avx_splitloop(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_flip_inline_modified(const int, const int, const int, const int, float *, float **, float *);

    //entropy search space reduction functions
    void find_optimal_mappings_scalar_replacement_modified_with_entropy(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_inline_modified_with_entropy(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_flip_ilp_modified_with_entropy(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_flip_ilp_modified_splitloop_with_entropy(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_avx_with_entropy(const int, const int, const int, const int, float *, float **, float *);
    void find_optimal_mappings_flip_avx_splitloop_with_entropy(const int, const int, const int, const int, float *, float **, float *);


    void find_optimal_mappings_flip_avx_all_transformations_in_vector(const int, const int, const int, const int, float *, float **, float *);
    
}


#endif // __FIC_FIND_OPTIMAL_MAPPINGS_H
