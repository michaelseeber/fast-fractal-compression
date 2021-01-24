cmake_minimum_required(VERSION 3.10)

# This file included below takes care of compiling dependencies
# and setting compiler flags on different platforms
include("CMakeConfig.txt")

# set the project name
project(fic)

include_directories(ext)

# add the executable
#add_executable(fic
#
#  # Header files
#  include/fic/common.h
#  include/fic/tsc_x86.h

#  # Source code files
#  src/main.cpp
#)

# base implementation = straight forward code
add_executable(fic_compress
  # Header files
  include/fic/common.h
  include/fic/rotate_vector.h
  include/fic/tsc_x86.h
  # Source code files
  src/block_store.cpp
  src/blur.cpp
  src/rotate_vector.cpp
  src/difference_norm.cpp
  src/map_one_matrix_to_another.cpp
  src/mapping_between_two_matrices.cpp
  src/compress.cpp
        src/blockwise_entropy.cpp)

#with flip
add_executable(fic_compress_flip
  # Header files
  include/fic/common.h
  include/fic/rotate_vector.h
  include/fic/hflip_vector.h
  include/fic/tsc_x86.h
  # Source code files
  src/block_store.cpp
  src/blur.cpp
  src/rotate_vector.cpp
  src/hflip_vector.cpp
  src/difference_norm.cpp
  src/map_one_matrix_to_another_flip.cpp
  src/mapping_between_two_matrices_flip.cpp
  src/compress_flip.cpp
        src/blockwise_entropy.cpp)

add_executable(fic_decompress
  # Header files
  include/fic/common.h
  # Source code files
  src/resize.cpp
  src/rotate.cpp
  src/split_into_blocks.cpp
  src/decompress.cpp
)

add_executable(benchmark_difference_norm
  # Header files
  include/fic/common.h
  include/fic/tsc_x86.h
  # Source code files
  src/difference_norm.cpp
  benchmark/benchmark_difference_norm.cpp
)

add_executable(benchmark_find_optimal_mappings
  # Header files
  include/fic/common.h
  include/fic/tsc_x86.h
  # Source code files
  src/difference_norm.cpp
  src/find_optimal_mappings.cpp
  src/blockwise_sum.cpp
  src/blockwise_sum_of_squares.cpp
  benchmark/benchmark_find_optimal_mappings.cpp
  src/blockwise_entropy.cpp)

add_executable(fic_compress_fast
  # Header files
  include/fic/common.h
  # Source code files
  src/resize.cpp
  src/rotate.cpp
  src/split_into_blocks.cpp
  src/difference_norm.cpp
  src/blockwise_sum.cpp
  src/blockwise_sum_of_squares.cpp
  src/blockwise_sum_of_xmuly.cpp
  src/find_optimal_mappings.cpp
  src/compress_fast.cpp
  src/blockwise_entropy.cpp)

add_executable(fic_compress_flip_fast_reference
  # Header files
  include/fic/common.h
  # Source code files
  src/resize.cpp
  src/rotate.cpp
  src/split_into_blocks.cpp
  src/difference_norm.cpp
  src/blockwise_sum.cpp
  src/blockwise_sum_of_squares.cpp
  src/find_optimal_mappings.cpp
  src/compress_flip_fast.cpp
  src/blockwise_entropy.cpp)
target_compile_definitions(fic_compress_flip_fast_reference PUBLIC _FAST_REFERENCE_)

add_executable(fic_compress_flip_fast_scalar_replacement
  # Header files
  include/fic/common.h
  # Source code files
  src/resize.cpp
  src/rotate.cpp
  src/split_into_blocks.cpp
  src/difference_norm.cpp
  src/blockwise_sum.cpp
  src/blockwise_sum_of_squares.cpp
  src/find_optimal_mappings.cpp
  src/compress_flip_fast.cpp
  src/blockwise_entropy.cpp)
target_compile_definitions(fic_compress_flip_fast_scalar_replacement PUBLIC _SCALAR_REPLACEMENT_)

add_executable(fic_compress_flip_fast_inline
  # Header files
  include/fic/common.h
  # Source code files
  src/resize.cpp
  src/rotate.cpp
  src/split_into_blocks.cpp
  src/difference_norm.cpp
  src/blockwise_sum.cpp
  src/blockwise_sum_of_squares.cpp
  src/find_optimal_mappings.cpp
  src/compress_flip_fast.cpp
  src/blockwise_entropy.cpp)
target_compile_definitions(fic_compress_flip_fast_inline PUBLIC _INLINE_)

add_executable(fic_compress_flip_fast_ilp
  # Header files
  include/fic/common.h
  # Source code files
  src/resize.cpp
  src/rotate.cpp
  src/split_into_blocks.cpp
  src/difference_norm.cpp
  src/blockwise_sum.cpp
  src/blockwise_sum_of_squares.cpp
  src/find_optimal_mappings.cpp
  src/compress_flip_fast.cpp
  src/blockwise_entropy.cpp)
target_compile_definitions(fic_compress_flip_fast_ilp PUBLIC _ILP_)

add_executable(fic_compress_flip_fast_avx
  # Header files
  include/fic/common.h
  # Source code files
  src/resize.cpp
  src/rotate.cpp
  src/split_into_blocks.cpp
  src/difference_norm.cpp
  src/blockwise_sum.cpp
  src/blockwise_sum_of_squares.cpp
  src/find_optimal_mappings.cpp
  src/compress_flip_fast.cpp
  src/blockwise_entropy.cpp)
target_compile_definitions(fic_compress_flip_fast_avx PUBLIC _MANUAL_AVX_)

add_executable(fic_compress_flip_fast_entropy
  # Header files
  include/fic/common.h
  # Source code files
  src/resize.cpp
  src/rotate.cpp
  src/split_into_blocks.cpp
  src/difference_norm.cpp
  src/blockwise_sum.cpp
  src/blockwise_sum_of_squares.cpp
  src/find_optimal_mappings.cpp
  src/compress_flip_fast.cpp
  src/blockwise_entropy.cpp)
target_compile_definitions(fic_compress_flip_fast_entropy PUBLIC _ENTROPY_)

add_executable(fic_compress_flip_fast
  # Header files
  include/fic/common.h
  # Source code files
  src/resize.cpp
  src/rotate.cpp
  src/split_into_blocks.cpp
  src/difference_norm.cpp
  src/blockwise_sum.cpp
  src/blockwise_sum_of_squares.cpp
  src/find_optimal_mappings.cpp
  src/compress_flip_fast.cpp
  src/blockwise_entropy.cpp)

#add_executable(stb_image_demo.cpp
#  # Header files
#  include/fic/common.h
#  # Source code files
#  src/stb_image_demo.cpp
#)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# specify CMAKE_BUILD_TYPE and Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math -funroll-loops")

# Link to dependency libraries
target_link_libraries(fic_compress ${extra_libs})

