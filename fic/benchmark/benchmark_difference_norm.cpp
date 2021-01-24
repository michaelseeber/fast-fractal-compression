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
#include "fic/difference_norm.h"

using namespace std;

#define NR 32
#define CYCLES_REQUIRED 1e8
#define REP 50
#define EPS (1e-3)

const int blockDimension = 8;

/* prototype of the function you need to optimize */
typedef void (*comp_func)(const int, float *, float *, float *);

//headers
void register_functions();
double get_perf_score(comp_func f);
double perf_test(comp_func f, string desc, int flops);

void add_function(comp_func f, string name, int flop);

/* Global vars, used to keep track of student functions */
vector<comp_func> userFuncs;
vector<string> funcNames;
vector<int> funcFlops;
int numFuncs = 0;

void rands(float *m, size_t row, size_t col)
{
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<float> dist(-1.0, 1.0);
    for (size_t i = 0; i < row * col; ++i)
        m[i] = dist(gen);
}

void build(float **a, int m, int n)
{
    *a = static_cast<float *>(aligned_alloc(32, m * n * sizeof(float))); //TODO: align this to 32
    rands(*a, m, n);
}

void destroy(float *m)
{
    free(m);
}

/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions()
{
    add_function(&(reference::difference_norm), "difference_norm_reference", 1);
    //add_function(&(fast::difference_norm), "difference_norm_fast0", 1);
    add_function(&(fast::difference_norm_scalar_replacement), "difference_norm_scalar_replacement", 1);
    add_function(&(fast::difference_norm_vectorized), "difference_norm_vectorized", 1);
    add_function(&(fast::difference_norm_vectorized_without_normalization), "difference_norm_vectorized_without_sqrt", 1);
    //add_function(&(fast::difference_norm_vectorized), "difference_norm_modified", 1);
}

double nrm_sqr_diff(float *x, float *y, int n)
{
    float nrm_sqr = 0.0;
    for (int i = 0; i < n; i++)
    {
        nrm_sqr += (x[i] - y[i]) * (x[i] - y[i]);
    }

    if (isnan(nrm_sqr))
    {
        nrm_sqr = INFINITY;
    }

    return nrm_sqr;
}

/*
* Registers a user function to be tested by the driver program. Registers a
* string description of the function as well
*/
void add_function(comp_func f, string name, int flops)
{
    userFuncs.push_back(f);
    funcNames.emplace_back(name);
    funcFlops.push_back(flops);

    numFuncs++;
}

/*
* Checks the given function for validity. If valid, then computes and
* reports and returns the number of cycles required per iteration
*/
double perf_test(comp_func f, string desc, int flops)
{
    double cycles = 0.;
    size_t num_runs = 100;
    double multiplier = 1;
    myInt64 start, end;

    float *A, *B;
    float sum;

    build(&A, blockDimension, blockDimension);
    build(&B, blockDimension, blockDimension);

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do
    {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++)
        {
            f(blockDimension, A, B, &sum);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);

    list<double> cyclesList;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++)
    {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i)
        {
            f(blockDimension, A, B, &sum);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cyclesList.push_back(cycles);
    }
    total_cycles /= REP;

    destroy(A);
    destroy(B);
    cyclesList.sort();
    cycles = total_cycles; //cyclesList.front();
    return round((100.0 * flops) / cycles) / 100.0;
}

int flops(int blockDimension, string funcName){
    int f;
    if(funcName == "difference_norm"){
        f = 3 * blockDimension * blockDimension + 8 + 5;
    } else if(funcName == "difference_norm_scalar_replacement"){
        f = 3 * blockDimension * blockDimension + 8 + 5;

    } else if(funcName == "difference_norm_scalar_replacement_modified"){
        f = 3 * blockDimension * blockDimension;

    }else if(funcName == "difference_norm_vectorized") {
        f = 3 * blockDimension * blockDimension + 7 + 8 + 5;

    }else if(funcName == "difference_norm_vectorized_without_normalization") {
        f = 3 * blockDimension * blockDimension + 4;

    }else{
        f = 3 * blockDimension * blockDimension + 3;
    }
    return f;
}


int memory_accesses(int blockDimension, string funcName){
    int f;
    if(funcName == "difference_norm"){
        f = 2 * blockDimension * blockDimension + 1;
    } else if(funcName == "difference_norm_scalar_replacement"){
        f = 2 * blockDimension * blockDimension + 1;

    } else if(funcName == "difference_norm_scalar_replacement_modified"){
        f = 1 + 2 * blockDimension * blockDimension;

    }else if(funcName == "difference_norm_vectorized") {
        f = 1 + 2 * blockDimension * blockDimension;

    }else if(funcName == "difference_norm_vectorized_without_normalization") {
        f = 1 + 2 * blockDimension * blockDimension;

    }else if(funcName == "difference_norm_vectorized_without_sqrt") {
        f = 1 + 2 * blockDimension * blockDimension;

    }else{
        f = 2 * blockDimension * blockDimension;
    }
    //adjust for the number of bytes per float
    return f * 4;
}

int main(int argc, char **argv)
{
    cout << "Starting program. ";
    double perf;
    int i;

    register_functions();

    if (numFuncs == 0)
    {
        cout << endl;
        cout << "No functions registered - nothing for driver to do" << endl;
        cout << "Register functions by calling register_func(f, name)" << endl;
        cout << "in register_funcs()" << endl;

        return 0;
    }
    cout << numFuncs << " functions registered." << endl;

    //Check validity of functions.
    float *A, *B;
    float sum, sum_old, sum_base;
    build(&A, blockDimension, blockDimension);
    build(&B, blockDimension, blockDimension);

    sum_old = sum;
    reference::difference_norm(blockDimension, A, B, &sum);
    sum_base = sum;

    for (i = 0; i < numFuncs; i++)
    {
        sum = sum_old;
        comp_func f = userFuncs[i];
        f(blockDimension, A, B, &sum);
        float error = nrm_sqr_diff(&sum, &sum_base, 1);

        if (error > EPS)
        {
            cout << error << endl;
            cout << "ERROR!!!!  the results for the " << i + 1 << "th function are different to the previous" << std::endl;
            // return 0;
        }
    }
    destroy(A);
    destroy(B);
    int fl = 0;
    int mem_ac = 0;
    for (i = 0; i < numFuncs; i++)
    {
        string name = funcNames[i];
        fl = flops(blockDimension, name);
        mem_ac = memory_accesses(blockDimension, name);
        perf = perf_test(userFuncs[i], name, fl);
        cout << endl
             << "Running: " << name << endl;
        cout << perf << " flops / cycles" << endl;
        cout << fl << " flops" <<endl;
        cout << mem_ac << " bytes"<<endl;

    }

    return 0;
}
