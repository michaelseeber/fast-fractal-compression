# ASL Project

## Report 

The report can be found [here](22_report.pdf).

## Directories
* `mannyray/`: this is the original repository from [github](https://github.com/mannyray/fractalMapping)
    * we have done one tiny change: bugfix rotation (used for both compression and decompression)
    * only used for the `decompress` function (well, not anymore)
    * this is NOT our reference implementation and NOT our base implementation
* `fic/`: our code
    * `src/`
        * `compress.cpp`: the base implementation = the straight forward C code
            * uses double-precision floating point numbers
        * `compress_fast.cpp`: uses our fastest implementation of all functions
            * uses single-precision floating point numbers
        * We compare the optimized versions (i.e. namespace `fast`) of the performance-critical functions to the functions of the `reference` namespace
        * the implementation in the `reference` namespace already uses the algorithmic optimizations (6x speedup)
            * this enables us to compare the critical function `find_optimal_mappings`
    * `benchmark/` validate and measure performance of the performance-critical functions
        * `find_optimal_mappings`
            * `difference_norm`
* `test/`: python scripts for plotting (and running the `fic` code)


## Compile:

    cd team022/fic/
    mkdir release
    cd release/
    cmake ..
    make
    
## Run one of the fully optimized compression executables:
    ./fic_compress_fast
    ./fic_compress_flip_fast
    ./fic_compress_flip_fast_entropy (adjustments of entropy threshold have to be done manually in src/find_optimal_mappings.cpp)
Example with custom image and stride=2 (stride has to be even and >= 2): 
    ./fic_compress_flip_fast ../images/256_cameraman.bmp 2 
    
## Run decompression executable
    ./fic_decompress <relative/path/to/compression/file> <name_of_output_png_image>
Works with fic compression that uses only 4 rotations as well as fic compression that uses 2 flips * 4 rotations = 8 transformations.

## Run validation & benchmark executables (used the code in main.cpp from the last homework) for each performance critical function
    ./benchmark_find_optimal_mappings
Extremely useful when we later start with applying the techinques from the lecture to make the code run fast.
Validates and benchmarks the following performance critical functions that are used in fic/src/compress.cpp:


## Plotting
Both scripts support logging, such that plots can be created without rerunning.
* for the report it's best to use the `.pdf` plot and import it as a figure in latex

### Benchmark (roofline)
* `test/benchmark.py` use this to generate roofline plots. 

#### plot from report
```bash
python3 test/benchmark.py
```
* say `n` when it asks you if you want to compile
* say `y` when it asks you if you want to create the report plot

#### plot from presentation
* Measurements have been done using random 512x512 images. They are stored in the following files
    * [test/results_benchmark/result_23.05_21_04.txt](test/results_benchmark/result_23.05_21_04.txt):
    all measurements, including noflip, flip and unused versions of `find_optimal_mappings`
    * [test/results_benchmark/result_23.05_21_05.txt](test/results_benchmark/result_23.05_21_05.txt):
    the 5 versions without flip that build on top of each other
    * [test/results_benchmark/result_23.05_21_06.txt](test/results_benchmark/result_23.05_21_06.txt):
    the 5 versions with flip that build on top of each other
* to execute the plot, use `python3 test/benchmark.py`,
then say `n` and then enter the `23.05_21_0x` and replace x by either 4, 5, or 6


### Performance
* `test/performance.py` use this to generate performance plots
*  You can specifiy the input sizes in code and the script supports noflip, flip and both.

* Measurements used in the report are stored in the following file:
    * [test/results_performance/result_merged.txt](test/results_performance/result_merged.txt)
    * to get these measurements, we used `test/performance_flip.py`

## performance critcal function: `find_optimal_mappings`
* fast unoptimized a.k.a. starting from scratch: `find_optimal_mappings`
* scalar replacement: `find_optimal_mappings_scalar_replacement_modified`
    * compute RSS instead of Euclidean distance
* inline: `find_optimal_mappings_inline_modified`
* ~~code motion~~ failed
* ~~strength reduction~~ failed
* ILP
    * without flip: `find_optimal_mappings_ilp_modified`
    * with flip: `find_optimal_mappings_flip_ilp_modified_splitloop`
* manual vectorization:
    * without flip: `find_optimal_mappings_avx`
    * with flip: `find_optimal_mappings_flip_avx_splitloop`
* ~~reordered memory layout~~ failed


## Euler setup (RIP)
### General
* make sure you are in ETH network (VPN)
* use `ssh -CY <nethz name>@euler.ethz.ch`
    * the `-Y` option is for opening display stuff like plots
* add the following lines to your `~/.bashrc` (on euler)
```bash
env2lmod # enables the "new stack" of modules
module load gcc/8.2.0 cmake/3.15.3 python/3.6.4
```
* don't forget to source it after you change it
```bash
source ~/.bashrc
```

### Python
* setup virtual environment for runnting tests
```bash
python3 -m venv test/virtenv
source test/virtenv/bin/activate
pip install -r test/requirements.txt
```
* make sure the python scripts use the correct system info file
```bash
cp test/system_info_euler.json test/system_info.json
```
* whenever you log in, you need to source the virtual environment again
```bash
source test/virtenv/bin/activate
```

### Mounting via `sshfs`
This enables you to have all the project files on your machine so that you can use whatever editor or IDE you prefer.

```bash
mkdir asl-project-on-euler
sshfs <nethz name>@euler.ethz.ch:/cluster/home/<nethz name>/asl-project asl-project-on-euler
```

Add the following line to your local `~/.bashrc` file to make your live easier:
This enablesz you to use the command `mount-asl-project` to mount the whole thing.
```bash
alias mount-asl-project='sshfs <nethz name>@euler.ethz.ch:/cluster/home/<nethz name>/asl-project asl-project-on-euler'
```
