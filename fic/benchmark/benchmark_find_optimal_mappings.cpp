/**
*      _________   _____________________  ____  ______
*     / ____/   | / ___/_  __/ ____/ __ \/ __ \/ ____/
*    / /_  / /| | \__ \ / / / /   / / / / / / / __/
*   / __/ / ___ |___/ // / / /___/ /_/ / /_/ / /___
*  /_/   /_/  |_/____//_/  \____/\____/_____/_____/
*
*  http://www.acl.inf.ethz.ch/teaching/fastcode
*  How to Write Fast Numerical Code 263-2300 - ETH Zurich
*  Copyright (C) 2019 
*                   Tyler Smith        (smitht@inf.ethz.ch) 
*                   Alen Stojanov      (astojanov@inf.ethz.ch)
*                   Gagandeep Singh    (gsingh@inf.ethz.ch)
*                   Markus Pueschel    (pueschel@inf.ethz.ch)
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program. If not, see http://www.gnu.org/licenses/.
*/
//#include "stdafx.h"

#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <random>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "fic/tsc_x86.h"
#include "fic/find_optimal_mappings.h"

using namespace std;

#define NR 32
#define CYCLES_REQUIRED 1e8
#define REP 50
#define EPS (1e-3)

// if set to true we use minimal memory access count
bool ONLY_COMPULRORY_MISSES = true;

// command line arguments
int height = 512;
int width = 512;
int num_rot = 4; // without flip: 4, with flip: 8
int stride = 16;

// define constants
const int blockDimension = 8;
const int mapStore = 5;
const int block_num_rows = 8;
const int block_num_cols = 8;

// inferred constants (depend on the ones defined above)
int block_num_pixels;
int img_gray_num_blocks_vertical;
int img_gray_num_blocks_horizontal;
int img_gray_num_blocks;
int height_small;
int width_small;
int stride_small;
int img_small_num_blocks_vertical;
int img_small_num_blocks_horizontal;
int img_small_num_blocks;

/* prototype of the function you need to optimize */
typedef void (*comp_func)(const int, const int, const int, const int, float *, float **, float *);

//headers
void register_functions();
double get_perf_score(comp_func f);
double perf_test(comp_func f, int flops);

void add_function(comp_func f, string name, int flop);

/* Global vars, used to keep track of student functions */
vector<comp_func> userFuncs;
vector<string> funcNames;
vector<int> funcFlops;
int numFuncs = 0;

void rands(float *m, size_t row, size_t col) {
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<float> dist(-1.0, 1.0);
    for (size_t i = 0; i < row * col; ++i)
        m[i] = dist(gen);
}

void build(float **a, int m, int n) {
    *a = static_cast<float *>(aligned_alloc(32, m * n * sizeof(float))); //TODO: align this to 32
    rands(*a, m, n);
}

void destroy(float *m) {
    free(m);
}

/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions() {
    add_function(&(reference::find_optimal_mappings), "find_optimal_mappings_reference", 1);
    //add_function(&(reference::find_optimal_mappings_with_entropy), "find_optimal_mappings_with_entropy_reference", 1);
    add_function(&(fast::find_optimal_mappings_scalar_replacement), "find_optimal_mappings_scalar_replacement", 1);
    add_function(&(fast::find_optimal_mappings_scalar_replacement_modified), "find_optimal_mappings_scalar_replacement_modified", 1);
    add_function(&(fast::find_optimal_mappings_inline_modified), "find_optimal_mappings_inline_modified", 1);
    add_function(&(fast::find_optimal_mappings_strength_reduction_modified), "find_optimal_mappings_strength_reduction_modified", 1);
    add_function(&(fast::find_optimal_mappings_code_motion_modified), "find_optimal_mappings_code_motion_modified", 1);
    add_function(&(fast::find_optimal_mappings_ilp_modified), "find_optimal_mappings_ilp_modified", 1);
    add_function(&(fast::find_optimal_mappings_ilp_modified_reordered), "find_optimal_mappings_ilp_modified_reordered", 1);
    add_function(&(fast::find_optimal_mappings_GCC_vectorization_modified), "find_optimal_mappings_GCC_vectorization_modified", 1);
#if __AVX2__
    add_function(&(fast::find_optimal_mappings_avx), "find_optimal_mappings_avx", 1);
    add_function(&(fast::find_optimal_mappings_avx_splitloop), "find_optimal_mappings_avx_splitloop", 1);
    add_function(&(fast::find_optimal_mappings_avx_hsuminlined), "find_optimal_mappings_avx_hsuminlined", 1);
    add_function(&(fast::find_optimal_mappings_avx_althsum), "find_optimal_mappings_avx_althsum", 1);
    add_function(&(fast::find_optimal_mappings_avx_reordered), "find_optimal_mappings_avx_reordered", 1);
#endif

    // what's this function?
    //add_function(&(fast::find_optimal_mappings_ilp_spatial_locality_32kb_8way), "find_optimal_mappings spatial locality", 1);

}
void register_functions_flip() {
    add_function(&(reference::find_optimal_mappings_flip), "find_optimal_mappings_flip_reference", 1);
    // add_function(&(reference::find_optimal_mappings_flip_with_entropy), "find_optimal_mappings_flip_with_entropy_reference", 1);
    add_function(&(fast::find_optimal_mappings_flip_scalar_replacement_modified), "find_optimal_mappings_flip_scalar_replacement_modified", 1);
    add_function(&(fast::find_optimal_mappings_flip_inline_modified), "find_optimal_mappings_flip_inline_modified", 1);
    add_function(&(fast::find_optimal_mappings_flip_ilp_modified), "find_optimal_mappings_flip_ilp_modified", 1);
    add_function(&(fast::find_optimal_mappings_flip_ilp_modified_reordered), "find_optimal_mappings_flip_ilp_modified_reordered", 1);
    add_function(&(fast::find_optimal_mappings_flip_ilp_modified_splitloop), "find_optimal_mappings_flip_ilp_modified_splitloop", 1);
    add_function(&(fast::find_optimal_mappings_flip_ilp_modified_splitloop2), "find_optimal_mappings_flip_ilp_modified_splitloop2", 1);
#if __AVX2__
    add_function(&(fast::find_optimal_mappings_flip_avx), "find_optimal_mappings_flip_avx", 1);
    add_function(&(fast::find_optimal_mappings_flip_avx_splitloop), "find_optimal_mappings_flip_avx_splitloop", 1);
    add_function(&(fast::find_optimal_mappings_flip_avx_all_transformations_in_vector), "find_optimal_mappings_flip_avx_all_transformations_in_vector", 1);
    
#endif
}

float nrm_sqr_diff(float *x, float *y, int n) {
    float nrm_sqr = 0.0;
    for (int i = 0; i < n; i++) {
        nrm_sqr += (x[i] - y[i]) * (x[i] - y[i]);
    }

    if (isnan(nrm_sqr)) {
        nrm_sqr = INFINITY;
    }

    return nrm_sqr;
}

// TODO: might need a nicer fix here
double validate_mappings(float *x, float *y, int n) {
    float diff_sum = 0.0;
    for (int i = 0; i < n; i+=mapStore)
    {
        diff_sum += std::abs(x[i + 1] - y[i + 1]);
        diff_sum += std::abs(x[i + 2] - y[i + 2]);
        diff_sum += std::abs(x[i + 3] - y[i + 3]);
        diff_sum += std::abs(x[i + 4] - y[i + 4]);
    }

    if (isnan(diff_sum)) {
        diff_sum = INFINITY;
    }
    return diff_sum;
}

/*
* Registers a user function to be tested by the driver program. Registers a
* string description of the function as well
*/
void add_function(comp_func f, string name, int flops) {
    userFuncs.push_back(f);
    funcNames.emplace_back(name);
    funcFlops.push_back(flops);

    numFuncs++;
}

/*
* Checks the given function for validity. If valid, then computes and
* reports and returns the number of cycles required per iteration
*/
double perf_test(comp_func f, long flops) {
    double cycles = 0.;
    size_t num_runs = 1;
    double multiplier = 1;
    myInt64 start, end;

    float **img_small_blocks_rotations = new float *[num_rot];
    for (int i = 0; i < num_rot; ++i) img_small_blocks_rotations[i] = new float[img_small_num_blocks * block_num_pixels];
    float *img_gray_blocks = new float[img_gray_num_blocks * block_num_pixels];
    float *mappings = new float[img_gray_num_blocks * mapStore];
    for(int i = 0; i < img_gray_num_blocks * mapStore; i++){
        mappings[i] = std::numeric_limits<float>::max();
    }

    for (int i = 0; i < num_rot; ++i) build(&img_small_blocks_rotations[i], img_small_num_blocks, block_num_pixels);
    build(&img_gray_blocks, img_gray_num_blocks, block_num_pixels);

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            f(block_num_rows, block_num_cols, img_gray_num_blocks,
                img_small_num_blocks, img_gray_blocks, img_small_blocks_rotations, mappings);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);

    list<double> cyclesList;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            f(block_num_rows, block_num_cols, img_gray_num_blocks,
                img_small_num_blocks, img_gray_blocks, img_small_blocks_rotations, mappings);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cyclesList.push_back(cycles);
    }
    total_cycles /= REP;

    for (int i = 0; i < num_rot; ++i) destroy(img_small_blocks_rotations[i]);
    destroy(img_gray_blocks);
    destroy(mappings);
    cyclesList.sort();
    cycles = total_cycles; //cyclesList.front();
    std::cout << "flops, cycles: " << flops << ", " << cycles << std::endl;
    return round((100.0 * flops) / cycles) / 100.0;
}

//Return the number of total bytes moved for specific functions
long long int memory_accesses(string func_name, int img_gray_num_blocks, int img_small_num_blocks, int rot, int block_num_pixels) {
    long long int acc = 1;
    if (ONLY_COMPULRORY_MISSES) {
        if (func_name == "find_optimal_mappings_reference" || func_name == "find_optimal_mappings_flip_reference" ||
                func_name == "find_optimal_mappings_scalar_replacement" ||
                func_name == "find_optimal_mappings_flip_scalar_replacement_modified" ||
                func_name == "find_optimal_mappings_scalar_replacement_modified" ||
                func_name == "find_optimal_mappings_inline_modified" ||
                func_name == "find_optimal_mappings_strength_reduction_modified" ||
                func_name == "find_optimal_mappings_code_motion_modified" ||
                func_name == "find_optimal_mappings_ilp_modified"||
                func_name == "find_optimal_mappings_GCC_vectorization_modified" ||
                func_name == "find_optimal_mappings_flip_scalar_replacement_modified"||
                func_name == "find_optimal_mappings_flip_inline_modified" ||
                func_name == "find_optimal_mappings_flip_ilp_modified" ||
                func_name == "find_optimal_mappings_ilp_modified_reordered" ||
                func_name == "find_optimal_mappings_flip_ilp_modified_reordered" ||
                func_name == "find_optimal_mappings_flip_ilp_modified_splitloop"||
                func_name == "find_optimal_mappings_flip_ilp_modified_splitloop2" ||
                func_name == "find_optimal_mappings_flip_avx_all_transformations_in_vector" ||
                func_name == "find_optimal_mappings_avx" ||
                func_name == "find_optimal_mappings_avx_hsuminlined" ||
                func_name == "find_optimal_mappings_avx_splitloop" ||
                func_name == "find_optimal_mappings_avx_althsum" ||
                func_name == "find_optimal_mappings_avx_reordered") {
            acc = 0;
            acc = acc + img_gray_num_blocks * block_num_pixels; // img_gray_blocks
            acc = acc + img_small_num_blocks * block_num_pixels * num_rot; // img_small_blocks
            acc = acc + img_gray_num_blocks * mapStore; // mappings
        }else if (func_name == "find_optimal_mappings_with_entropy" ) {
            // I think this is not used
            acc = //elements accessed in every iteration over img_small_num_blocks
                    //small_blockwise_sumofsquares and small_blockwise_sums
                    ((rot + 1) * block_num_pixels * img_small_num_blocks +
                    //elements accessed in every iteration over img_gray_num_blocks
                    // img_small_blocks_rotations, img_gray_blocks, vec1, and mappings
                    + block_num_pixels * 2 + 5) * img_gray_num_blocks;
            //blockwise_sum memory accesses

            acc += img_gray_num_blocks * block_num_pixels * 2;
            acc += img_small_num_blocks * (2 + block_num_pixels);
            acc += img_small_num_blocks * (block_num_pixels + 1);
        } else if (func_name == "find_optimal_mappings_flip_avx" ||
            func_name == "find_optimal_mappings_flip_avx_splitloop") {
            // TODO: check
            acc = 0;
            acc = acc + img_gray_num_blocks * block_num_pixels; // img_gray_blocks
            acc = acc + img_small_num_blocks * block_num_pixels * num_rot; // img_small_blocks
            acc = acc + img_gray_num_blocks * mapStore; // mappings
        } else {
            cout << "function name not recognized number of memory acesses. using same acc as the reference function" << endl;
            //elements accessed in every iteration over img_small_num_blocks
                    //small_blockwise_sumofsquares and small_blockwise_sums
            acc = rot * block_num_pixels + 2;
            acc = acc * img_small_num_blocks;
                    //elements accessed in every iteration over img_gray_num_blocks
                    // img_small_blocks_rotations, img_gray_blocks, vec1, and mappings
            acc = acc + block_num_pixels + 5;
            acc = acc * img_gray_num_blocks;
            //blockwise_sum memory accesses
            acc += img_gray_num_blocks * block_num_pixels * 2;
            acc += img_small_num_blocks * (2 + block_num_pixels);
        }
    } else {
        if (func_name == "find_optimal_mappings_reference" || func_name == "find_optimal_mappings_flip_reference" ||
                func_name == "find_optimal_mappings_scalar_replacement" ||
                func_name == "find_optimal_mappings_flip_scalar_replacement_modified" ||
                func_name == "find_optimal_mappings_scalar_replacement_modified" ||
                func_name == "find_optimal_mappings_inline_modified" ||
                func_name == "find_optimal_mappings_strength_reduction_modified" ||
                func_name == "find_optimal_mappings_code_motion_modified" ||
                func_name == "find_optimal_mappings_ilp_modified"||
                func_name == "find_optimal_mappings_GCC_vectorization_modified" ||
                func_name == "find_optimal_mappings_flip_scalar_replacement_modified"||
                func_name == "find_optimal_mappings_flip_inline_modified" ||
                func_name == "find_optimal_mappings_flip_ilp_modified" ||
                func_name == "find_optimal_mappings_ilp_modified_reordered" ||
                func_name == "find_optimal_mappings_flip_ilp_modified_reordered" ||
                func_name == "find_optimal_mappings_flip_ilp_modified_splitloop"||
                func_name == "find_optimal_mappings_flip_ilp_modified_splitloop2" ||
                func_name == "find_optimal_mappings_flip_avx_all_transformations_in_vector") {
            //elements accessed in every iteration over img_small_num_blocks
                    //small_blockwise_sumofsquares and small_blockwise_sums
            acc = rot * block_num_pixels + 2;
            acc = acc * img_small_num_blocks;
                     //elements accessed in every iteration over img_gray_num_blocks
                     // img_small_blocks_rotations, img_gray_blocks, vec1, and mappings
            acc = acc + block_num_pixels + 5;
            acc = acc * img_gray_num_blocks;
            //blockwise_sum memory accesses
            acc += img_gray_num_blocks * block_num_pixels * 2;
            acc += img_small_num_blocks * (2 + block_num_pixels);
        }else if (func_name == "find_optimal_mappings_with_entropy" ) {
            acc = //elements accessed in every iteration over img_small_num_blocks
                    //small_blockwise_sumofsquares and small_blockwise_sums
                    ((rot + 1) * block_num_pixels * img_small_num_blocks +
                     //elements accessed in every iteration over img_gray_num_blocks
                     // img_small_blocks_rotations, img_gray_blocks, vec1, and mappings
                     + block_num_pixels * 2 + 5) * img_gray_num_blocks;
            //blockwise_sum memory accesses

            acc += img_gray_num_blocks * block_num_pixels * 2;
            acc += img_small_num_blocks * (2 + block_num_pixels);
            acc += img_small_num_blocks * (block_num_pixels + 1);
        } else if (func_name == "find_optimal_mappings_avx" ||
            func_name == "find_optimal_mappings_avx_hsuminlined" ||
            func_name == "find_optimal_mappings_avx_splitloop" ||
            func_name == "find_optimal_mappings_avx_althsum") {
            acc = img_small_num_blocks * (2 + 8 * 4 * 8 );
            acc =  img_gray_num_blocks * (1 + acc + 5 + 1);
            acc +=   //blockwise_sum memory accesses
                    img_gray_num_blocks * (block_num_pixels + 1)+
                     (img_small_num_blocks * (block_num_pixels + 2));
        } else if (func_name == "find_optimal_mappings_flip_avx" ||
            func_name == "find_optimal_mappings_flip_avx_splitloop") {
            // TODO: check
            acc = img_small_num_blocks * (2 + 8 * 4 * 8 );
            acc =  img_gray_num_blocks * (1 + acc + 5 + 1);
            acc +=   //blockwise_sum memory accesses
                    img_gray_num_blocks * (block_num_pixels + 1)+
                     (img_small_num_blocks * (block_num_pixels + 2));
        } else if (func_name == "find_optimal_mappings_avx_reordered") {
            acc = img_small_num_blocks * (2 + 8 * 5 * 8 +
                                          8 * 5 * 8 + 5);
            acc =   //blockwise_sum memory accesses
                    img_gray_num_blocks * block_num_pixels * 2 +
                    2 * (img_small_num_blocks * block_num_pixels * 2) +
                    img_gray_num_blocks * (1 + acc);
            acc += img_small_num_blocks * 8 * 3 * rot * 2;
        } else {
            cout << "function name not recognized number of memory acesses. using same acc as the reference function" << endl;
            //elements accessed in every iteration over img_small_num_blocks
                    //small_blockwise_sumofsquares and small_blockwise_sums
            acc = rot * block_num_pixels + 2;
            acc = acc * img_small_num_blocks;
                     //elements accessed in every iteration over img_gray_num_blocks
                     // img_small_blocks_rotations, img_gray_blocks, vec1, and mappings
            acc = acc + block_num_pixels + 5;
            acc = acc * img_gray_num_blocks;
            //blockwise_sum memory accesses
            acc += img_gray_num_blocks * block_num_pixels * 2;
            acc += img_small_num_blocks * (2 + block_num_pixels);
        }
    }
    // adjust for the number of bytes per float
    return acc * 4;
}

    void print_usage() {
        std::cout << "usage: ./benchmark_find_optimal_mappings [<width>[x<height>] [<'flip' or 'noflip'> [<stride>]]]"
                  << endl;
        std::cout << "examples:" << endl;
        std::cout << "\t./benchmark_find_optimal_mappings" << endl;
        std::cout << "\t./benchmark_find_optimal_mappings 256" << endl;
        std::cout << "\t./benchmark_find_optimal_mappings 256x512" << endl;
        std::cout << "\t./benchmark_find_optimal_mappings 256 flip" << endl;
        std::cout << "\t./benchmark_find_optimal_mappings 256x512 noflip 8" << endl;
        std::cout << "The default is:" << endl;
        std::cout << "\t./benchmark_find_optimal_mappings 512x512 noflip 16" << endl;
    }

    int main(int argc, char **argv) {

        // first command line argument: image size
        if (argc > 1) {
            // check if width and height is given or only one of them
            char *c = argv[1];
            char *width_in = new char[8];
            char *height_in = new char[8];
            int index = 0;
            while (*c) {
                if (strchr("x", *c)) {
                    // both width and height are given
                    width = atoi(width_in);
                    c++;
                    break;
                } else {
                    width_in[index] = *c;
                }
                c++;
                index++;
            }
            if (*c) {
                index = 0;
                while (*c) {
                    height_in[index] = *c;
                    index++;
                    c++;
                }
                height = atoi(height_in);
            } else {
                width = atoi(width_in);
                height = width;
            }
            if (height % 8 || width % 8) {
                std::cout << "invalid resoltion: " << width << "x" << height << ". please use multiples of 8." << endl;
                print_usage();
                std::cout << "exiting" << endl;
                exit(-1);
            }
        }

        // second command line argument: flip or noflip
        if (argc > 2) {
            if (!strcmp(argv[2], "flip") || !strcmp(argv[2], "noflip")) {
                num_rot = !strcmp(argv[2], "flip") ? 8 : 4;
            } else {
                std::cout << "second command line argument must be either \"flip\" or \"noflip\"" << endl;
                print_usage();
                std::cout << "exiting" << endl;
                exit(-1);
            }
        }

        // third command line argument: stride
        if (argc > 3) {
            stride = atoi(argv[3]);
            if (stride != 2 && stride != 4 && stride != 8 && stride != 16) {
                std::cout << "stride must be one of {2, 4, 8, 16}" << endl;
                print_usage();
                std::cout << "exiting" << endl;
                exit(-1);
            }
        } // done reading command line arguments

        std::cout << "Starting program for random image of size " << width << " x " << height << std::endl;;

        block_num_pixels = block_num_cols * block_num_rows;
        img_gray_num_blocks_vertical = height / block_num_cols;
        img_gray_num_blocks_horizontal = width / block_num_rows;
        img_gray_num_blocks = img_gray_num_blocks_vertical * img_gray_num_blocks_horizontal;
        height_small = height / 2;
        width_small = width / 2;
        stride_small = stride / 2;
        img_small_num_blocks_vertical = height_small / stride_small - (block_num_cols / stride_small) + 1;
        img_small_num_blocks_horizontal = width_small / stride_small - (block_num_rows / stride_small) + 1;
        img_small_num_blocks = img_small_num_blocks_vertical * img_small_num_blocks_horizontal;
        double perf;
        int i;

        if (num_rot == 4) {
            // without flip
            register_functions();
        } else {
            // with flip
            register_functions_flip();
        }

        if (numFuncs == 0) {
            std::cout << endl;
            std::cout << "No functions registered - nothing for driver to do" << endl;
            std::cout << "Register functions by calling register_func(f, name)" << endl;
            std::cout << "in register_funcs()" << endl;

            return 0;
        }
        std::cout << numFuncs << " functions registered." << endl;

        //Check validity of functions.
        float **img_small_blocks_rotations = new float *[num_rot];
        for (int i = 0; i < num_rot; ++i)
            img_small_blocks_rotations[i] = new float[img_small_num_blocks * block_num_pixels];
        float *img_gray_blocks = new float[img_gray_num_blocks * block_num_pixels];
        float *mappings = new float[img_gray_num_blocks * mapStore];

        for (int i = 0; i < num_rot; ++i) build(&img_small_blocks_rotations[i], img_small_num_blocks, block_num_pixels);
        build(&img_gray_blocks, img_gray_num_blocks, block_num_pixels);

        float *mappings_reference = new float[img_gray_num_blocks * mapStore];
        for (int i = 0; i < img_gray_num_blocks * mapStore; i++) {
            mappings_reference[i] = std::numeric_limits<float>::max();
        }
        if (num_rot == 4) {
            // without flip
            reference::find_optimal_mappings(block_num_rows, block_num_cols,
                                             img_gray_num_blocks,
                                             img_small_num_blocks, img_gray_blocks, img_small_blocks_rotations,
                                             mappings_reference);
        } else {
            // with flip
            reference::find_optimal_mappings_flip(block_num_rows, block_num_cols,
                                                  img_gray_num_blocks,
                                                  img_small_num_blocks, img_gray_blocks, img_small_blocks_rotations,
                                                  mappings_reference);
        }

        for (i = 0; i < numFuncs; i++) {
            for (int i = 0; i < img_gray_num_blocks * mapStore; i++) {
                mappings[i] = std::numeric_limits<float>::max();
            }
            comp_func f = userFuncs[i];
            f(block_num_rows, block_num_cols, img_gray_num_blocks,
              img_small_num_blocks, img_gray_blocks, img_small_blocks_rotations, mappings);
            double error = validate_mappings(mappings, mappings_reference, img_gray_num_blocks * mapStore);

            if (error > EPS) {
                std::cout << "ERROR!!!!  the results for function '" << funcNames[i]
                          << "' are different to the previous" << std::endl;
                std::cout << "error = " << error << endl;
                // return 0;
            }
        }
        for (int i = 0; i < num_rot; ++i) destroy(img_small_blocks_rotations[i]);
        destroy(img_gray_blocks);
        destroy(mappings);

        // compute number of flops (incl divs and sqrts)
        long long int num_flops = 0;
        long num_divs = 0;
        long num_sqrts = 0;
        for (i = 0; i < numFuncs; i++) {
            num_flops = 1;
            num_divs = 1;
            num_sqrts = 1;
            if (funcNames[i] == "find_optimal_mappings_reference" ||
                funcNames[i] == "find_optimal_mappings_flip_reference" ||
                funcNames[i] == "find_optimal_mappings_scalar_replacement") {

                // nested loop
                num_flops = num_flops * 2 * block_num_pixels; // compute graymulsmall_blockwise_sum
                num_flops = num_flops + 3 + 3; // compute contrast & brightness
                num_flops = num_flops + 2 * block_num_pixels; // compute vec1
                num_flops = num_flops + block_num_pixels * 3; // difference norm
                num_flops = num_flops * num_rot; // 3rd layer loop
                num_flops = num_flops + 3; // compute denom
                num_flops = num_flops * img_small_num_blocks; // 2nd layer loop
                num_flops = num_flops * img_gray_num_blocks; // 1st layer loop

                // precomputations
                num_flops = num_flops + img_gray_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels * 2;

                // number of divisions
                num_divs = 1; // contrast
                num_divs = num_divs + 1; // brightness
                num_divs = num_divs + 2; // difference_norm
                num_divs = num_divs * num_rot * img_small_num_blocks * img_gray_num_blocks; // inner 3 loops

                // number of square roots
                num_sqrts = 1; // the one sqrt in difference norm
                num_sqrts = num_sqrts * num_rot * img_small_num_blocks * img_gray_num_blocks; // inner 3 loops

                std::cout << num_flops << " flops + " << num_divs << " divs + " << num_sqrts << " sqrts" << endl;
            } else if (funcNames[i] == "find_optimal_mappings_scalar_replacement_modified" ||
                funcNames[i] == "find_optimal_mappings_flip_scalar_replacement_modified") {
                // nested loop
                num_flops = num_flops * 2 * block_num_pixels; // compute graymulsmall_blockwise_sum
                num_flops = num_flops + 3 + 3; // compute contrast & brightness
                num_flops = num_flops + 2 * block_num_pixels; // compute vec1
                num_flops = num_flops + block_num_pixels * 3; // difference norm
                num_flops = num_flops * num_rot; // 3rd layer loop
                num_flops = num_flops + 3; // compute denom
                num_flops = num_flops * img_small_num_blocks; // 2nd layer loop
                num_flops = num_flops * img_gray_num_blocks; // 1st layer loop

                // precomputations
                num_flops = num_flops + img_gray_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels * 2;

                // number of divisions
                num_divs = 1; // contrast
                num_divs = num_divs + 1; // brightness
                num_divs = num_divs * num_rot * img_small_num_blocks * img_gray_num_blocks; // inner 3 loops

                // number of square roots
                num_sqrts = 0;
            } else if (funcNames[i] == "find_optimal_mappings_inline_modified" ||
                funcNames[i] == "find_optimal_mappings_flip_inline_modified") {
                
                // nested loop
                num_flops = num_flops * 2 * block_num_pixels; // compute graymulsmall_blockwise_sum
                num_flops = num_flops + 3 + 3; // compute contrast & brightness
                num_flops = num_flops + block_num_pixels * (2 + 1 + 1 + 1);
                num_flops = num_flops * num_rot; // 3rd layer loop
                num_flops = num_flops + 3; // compute denom
                num_flops = num_flops * img_small_num_blocks; // 2nd layer loop
                num_flops = num_flops * img_gray_num_blocks; // 1st layer loop

                // precomputations
                num_flops = num_flops + img_gray_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels * 2;

                // number of divisions
                num_divs = 1; // contrast
                num_divs = num_divs + 1; // brightness
                num_divs = num_divs * num_rot * img_small_num_blocks * img_gray_num_blocks; // inner 3 loops

                // number of square roots
                num_sqrts = 0;
            } else if (funcNames[i] == "find_optimal_mappings_strength_reduction_modified") {
                // NOT an optimization
                // nested loop
                num_flops = num_flops * 2 * block_num_pixels; // compute graymulsmall_blockwise_sum
                num_flops = num_flops + 4 + 4; // compute contrast & brightness (now with a multiplication)
                num_flops = num_flops + block_num_pixels * (2 + 1 + 1 + 1);
                num_flops = num_flops * num_rot; // 3rd layer loop
                num_flops = num_flops + 3;
                num_flops = num_flops * img_small_num_blocks; // 2nd layer loop
                num_flops = num_flops * img_gray_num_blocks; // 1st layer loop

                // precomputations
                num_flops = num_flops + img_gray_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels * 2;

                // number of divisions
                num_divs = 1; // one division before looping over rotations
                num_divs = num_divs * img_small_num_blocks * img_gray_num_blocks; // inner 2 loops

                // number of square roots
                num_sqrts = 0;
            } else if (funcNames[i] == "find_optimal_mappings_code_motion_modified") {
                // NOT an optimization
                // nested loop
                num_flops = num_flops * 2 * block_num_pixels; // compute graymulsmall_blockwise_sum
                num_flops = num_flops + 2 + 2; // compute contrast & brightness (with one flop less each)
                num_flops = num_flops + block_num_pixels * (2 + 1 + 1 + 1);
                num_flops = num_flops * num_rot; // 3rd layer loop
                num_flops = num_flops + 3; // compute denom
                num_flops = num_flops + 2; // compute temps (mult only)
                num_flops = num_flops * img_small_num_blocks; // 2nd layer loop
                num_flops = num_flops * img_gray_num_blocks; // 1st layer loop

                // precomputations
                num_flops = num_flops + img_gray_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels * 2;

                // number of divisions
                num_divs = 1; // temp0
                num_divs = num_divs + 1; // temp1
                num_divs = num_divs + 1; // temp2
                num_divs = num_divs * img_small_num_blocks * img_gray_num_blocks; // inner 3 loops

                // number of square roots
                num_sqrts = 0;

            } else if (funcNames[i] == "find_optimal_mappings_ilp_modified" ||
                funcNames[i] == "find_optimal_mappings_ilp_modified_reordered" ||
                funcNames[i] == "find_optimal_mappings_flip_ilp_modified" ||
                funcNames[i] == "find_optimal_mappings_flip_ilp_modified_splitloop" ||
                funcNames[i] == "find_optimal_mappings_flip_ilp_modified_splitloop2" ||
                funcNames[i] == "find_optimal_mappings_flip_ilp_modified_reordered") {
                
                // nested loop
                num_flops = 3; // compute denom
                num_flops = num_flops + num_rot * 2 * block_num_pixels; // compute graymulsmall_blockwise_sum
                num_flops = num_flops + num_rot * 3; // compute contrast
                num_flops = num_flops + num_rot * 3; // compute brightness
                num_flops = num_flops + block_num_pixels * num_rot * (2 + 1 + 1 + 1);
                num_flops = num_flops * img_small_num_blocks; // 2nd layer loop
                num_flops = num_flops * img_gray_num_blocks; // 1st layer loop

                // precomputations
                num_flops = num_flops + img_gray_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels * 2;

                // number of divisions
                num_divs = num_rot; // contrast
                num_divs = num_divs + num_rot; // brightness
                num_divs = num_divs * img_small_num_blocks * img_gray_num_blocks; // inner 2 loops

                // number of square roots
                num_sqrts = 0;
            } else if (funcNames[i] == "find_optimal_mappings_avx" ||
                funcNames[i] == "find_optimal_mappings_avx_splitloop" ||
                funcNames[i] == "find_optimal_mappings_avx_hsuminlined" ||
                funcNames[i] == "find_optimal_mappings_flip_avx" ||
                funcNames[i] == "find_optimal_mappings_flip_avx_splitloop" ||
                funcNames[i] == "find_optimal_mappings_avx_althsum" ||
                funcNames[i] == "find_optimal_mappings_avx_reordered") {

                
                // nested loop
                num_flops = 3; // compute denom
                num_flops = num_flops + 8 * 8 * num_rot * 2; // compute graymulsmall...
                //num_flops = num_flops + 3 * 8 * num_rot; // compute horizontal_sum
                num_flops = num_flops + 7 * num_rot; // compute horizontal_sum
                num_flops = num_flops + 3 * num_rot; // compute contrast
                num_flops = num_flops + 3 * num_rot; // compute brightness
                num_flops = num_flops + 8 * (8 * 2 * num_rot + 8 * num_rot + 8 * 2 * num_rot); // compute error
                //num_flops = num_flops + 3 * 8 * num_rot; // compute horizontal_sum again
                num_flops = num_flops + 7 * num_rot; // compute horizontal_sum again
                num_flops = num_flops * img_small_num_blocks; // 2nd layer loop
                num_flops = num_flops * img_gray_num_blocks; // 1st layer loop

                // inlined blockwise sum
                num_flops = num_flops + img_gray_num_blocks * block_num_pixels; // compute gray_blockwise_sums
                num_flops = num_flops + 3 * img_small_num_blocks * block_num_pixels; // compute small_blockwise_sums and small_blockwise_sumofsquares

                // number of divisions
                num_divs = num_rot; // contrast
                num_divs = num_divs + num_rot; // brightness
                num_divs = num_divs * img_small_num_blocks * img_gray_num_blocks; // inner 2 loops

                // number of square roots
                num_sqrts = 0;
            } else if (funcNames[i] == "find_optimal_mappings_flip_avx_all_transformations_in_vector")
            {
                // nested loop
                num_flops = 3; // compute denom
                num_flops = num_flops + block_num_pixels * 2 *8; //compute graymulsmall
                num_flops = num_flops + 8 * 4; //fma
                num_flops = num_flops + 8 * 2; //mults
                num_flops = num_flops + block_num_pixels * 5 *8; //compute error
                num_flops = num_flops + 4; //set1 precompute
                num_flops = num_flops * img_small_num_blocks; // 2nd layer loop
                num_flops = num_flops * img_gray_num_blocks; // 1st layer loop

                // inlined blockwise sum
                num_flops = num_flops + img_gray_num_blocks * block_num_pixels; // compute gray_blockwise_sums
                num_flops = num_flops + 3 * img_small_num_blocks * block_num_pixels; // compute small_blockwise_sums and small_blockwise_sumofsquares

                num_divs = 1*img_small_num_blocks*img_gray_num_blocks; //denom
            } else if (funcNames[i] == "find_optimal_mappings_GCC_vectorization_modified") {
                
                // nested loop
                num_flops = num_flops * 2 * block_num_pixels; // compute graymulsmall_blockwise_sum
                num_flops = num_flops + 3 + 3; // compute contrast & brightness
                num_flops = num_flops + block_num_pixels * (2 + 1 + 1 + 1);
                num_flops = num_flops * num_rot; // 3rd layer loop
                num_flops = num_flops + 3; // compute denom
                num_flops = num_flops * img_small_num_blocks; // 2nd layer loop
                num_flops = num_flops * img_gray_num_blocks; // 1st layer loop

                // precomputations
                num_flops = num_flops + img_gray_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels;
                num_flops = num_flops + img_small_num_blocks * block_num_pixels * 2;

                // number of divisions
                num_divs = 1; // contrast
                num_divs = num_divs + 1; // brightness
                num_divs = num_divs * num_rot * img_small_num_blocks * img_gray_num_blocks; // inner 3 loops

                // number of square roots
                num_sqrts = 0;
            } else {
                std::cout << "function name not recognized for flop count. using 0 flops" << endl;
            }
            long total_num_flops = num_flops + num_divs + num_sqrts;
            long long int total_num_bytes = memory_accesses(funcNames[i], img_gray_num_blocks, img_small_num_blocks, num_rot,
                                                            block_num_pixels);

            perf = perf_test(userFuncs[i], total_num_flops);
            cout << "Running: " << funcNames[i] << endl;
            cout << perf << " flops / cycles" << endl;
            cout << total_num_flops << " flops" << endl;
            cout << total_num_bytes << " bytes" << endl;
            if (num_divs) {
                std::cout << ((float) num_divs) / total_num_flops * 100. << "% of them are divisions" << endl;
            }
            if (num_sqrts) {
                std::cout << ((float) num_sqrts) / total_num_flops * 100. << "% of them are square roots" << endl;
            }

            std::cout << endl;
        }

        return 0;
    }

