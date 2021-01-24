import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import time
import re
import json

from subprocess import Popen, PIPE


def run_benchmark(benchmark):
    p = Popen(benchmark["cmd"], cwd=working_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    print("Running Benchmark " + benchmark["name"] +" ...")
    benchmark_output, benchmark_errors = p.communicate()

    if benchmark_errors != '':
        print(benchmark_errors)
    else:
        return benchmark_output


#Constants and arguments
try:
    with open('test/system_info.json') as f:
        system_info_dict = json.load(f)
except Exception as e:
    print(e)
    print('')
    print('Please create your system info file. You can use test/system_info_michi.json as an example. Then copy it to "test/system_info.json"')
    exit(1)
SYSTEM_INFO = "" + system_info_dict['cpu_name'] + ' @ ' + system_info_dict['cpu_base_freq'] \
    + '\nCompiler: ' + system_info_dict['compiler'] + ', OS: ' + system_info_dict['os']

difference_norm = { "name" : "difference_norm", "cmd" : "./benchmark_difference_norm"}
find_optimal_mapping= { "name": "optimal_mappings", "cmd" : ["./benchmark_find_optimal_mappings", "1024", "noflip"]}
compile_cmd = "make"
working_dir = 'fic/release/'

#load from log or recompile
answer = input("Yes to compile and run benchmarks. No to plot form logs. (y/n) ")
print("that was the answer")
print(answer)
if answer == "y":
    #compile code
    print("Compiling...")
    p = Popen(compile_cmd, cwd=working_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    compile_output, compile_errors = p.communicate()

    if compile_errors != '':
        print(compile_errors)
        answer = input("Warnings or Errors detected, do you still want to continue? (y/n) ")
        if answer == "n":
            exit(1)
    else:
        print(compile_output)

    
    image_size = input("Enter size to use for find_optimal_mappings (256, 512, 1024): ")
    find_optimal_mapping["cmd"][1] = str(image_size)
    flip = input("Enter either 'flip' or 'noflip' or 'all': ")
    runall = False
    if(flip == "all"):
        runall = True
    else:
        find_optimal_mapping["cmd"][2] = str(flip)
    #benchmark code
    print("Benchmarking...")
    benchmarks_result = run_benchmark(difference_norm)
    benchmarks_result += run_benchmark(find_optimal_mapping)
    if(runall):
        find_optimal_mapping["cmd"][2] = "flip" #for all we also run with flip
        benchmarks_result += run_benchmark(find_optimal_mapping)
    #store for logging
    path_to_output_file = 'test/results_benchmark/result_' + time.strftime("%d.%m_%H:%M") +'.txt'
    with open(path_to_output_file, 'w+') as output_file:
        output_file.write(benchmarks_result)
elif answer == "n":

    report = input("You want to create the report plot?(y/n) ")

    if report == "y":
        path_to_log_file_all = 'test/results_benchmark/result_23.05_21_04.txt'
        path_to_log_file_noflip = 'test/results_benchmark/result_23.05_21_05.txt'
        path_to_log_file_flip = 'test/results_benchmark/result_23.05_21_06.txt'
        try:
            with open(path_to_log_file_all, 'r') as log:
                benchmarks_result = log.read()
            with open(path_to_log_file_noflip, 'r') as log:
                benchmarks_result_noflip = log.read()
            with open(path_to_log_file_flip, 'r') as log:
                benchmarks_result_flip = log.read()
        
        except IOError:
            print ("Could not read files to plot report plots.")
            exit(1)
            

    elif report == "n":
        log_date = input("Enter time of log. Example: 26.04_15:40  ")
        path_to_log_file = 'test/results_benchmark/result_' + log_date +'.txt'
        try:
            with open(path_to_log_file, 'r') as log:
                benchmarks_result = log.read()

        except IOError:
            print ("Could not read file: " + path_to_log_file)
            exit(1)

    else:
        print("No valid input")
        exit(1)

    
else:
    print("No valid input")
    exit(1)



#plotting code
parse = re.findall("Running:\s*(\S*)\s(\d*.?\d*)\s\S*\s\S\s\S*\s(\d*)\s\S*\s(\d*)", benchmarks_result)

#rooflineplot
sns.set()
# sns.set_palette(sns.hls_palette(8, l=.3, s=.8), color_codes=True)

fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(1, 1, 1)

x_start = 0.1
x_end = 1000000
#roofs
I_vector = [x_start, x_end]
ax.plot(I_vector, list(map(lambda x: system_info_dict['peak_perf_vec'], I_vector)), 'b-', color='forestgreen', label='Peak Performance \u03A0 (vector)')
ax.plot(I_vector, list(map(lambda x: system_info_dict['peak_perf'], I_vector)), 'b-', color='mediumseagreen', label='Peak Performance \u03A0 (scalar)')
ax.plot(I_vector, list(map(lambda x: system_info_dict['mem_bandwidth_theoretical']*x, I_vector)), 'b-', color='maroon', label='Read bandwidth theoretical \u03B2')
ax.plot(I_vector, list(map(lambda x: system_info_dict['mem_bandwidth_measured']*x, I_vector)), 'b-', color='indianred', label='Read bandwidth measured \u03B2')


for datapoint in parse: #add datapoints
    operational_intensity = float(datapoint[2]) / float(datapoint[3])
    print(operational_intensity)
    plt.scatter(operational_intensity,float(datapoint[1]),marker='o',color='lightblue')
    #plt.text(operational_intensity+0.005,float(datapoint[1]),s=index)


if benchmarks_result_noflip != None and benchmarks_result_flip != None:
    parsenoflip = re.findall("Running:\s*(\S*)\s(\d*.?\d*)\s\S*\s\S\s\S*\s(\d*)\s\S*\s(\d*)", benchmarks_result_noflip)
    parseflip = re.findall("Running:\s*(\S*)\s(\d*.?\d*)\s\S*\s\S\s\S*\s(\d*)\s\S*\s(\d*)", benchmarks_result_flip)

    for datapoint in parsenoflip: #
        operational_intensity = float(datapoint[2]) / float(datapoint[3])
        print(operational_intensity)
        plt.scatter(operational_intensity,float(datapoint[1]),marker='o',color='royalblue')
        #plt.text(operational_intensity+0.005,float(datapoint[1]),s=index)

    annotation_function_names =["Fast Unoptimized", "Scalar Replacement", "Inlining", "ILP", "Manual Vectorization" ]
    index = 0
    for datapoint in parseflip: #n
        operational_intensity = float(datapoint[2]) / float(datapoint[3])
        print(operational_intensity)
        plt.scatter(operational_intensity,float(datapoint[1]),marker='o',color='navy')


        # plt.text(operational_intensity+3000,float(datapoint[1]),s=annotation_function_names[index], transform=ax.transAxe)

        ax.annotate(annotation_function_names[index], xy=(operational_intensity,float(datapoint[1])), xytext=(10,-5), textcoords="offset points")
        index +=1

   


ax.set_title(SYSTEM_INFO, loc="left", pad=30, fontweight="bold")
ax.set_xlabel('Operational Intensity [Flops/Byte]')
ax.set_ylabel('Performance [Flops/Cycle]', rotation=0)
ax.yaxis.set_label_coords(0.15, 1.02)

ax.set_ylim(0.6, 64)
ax.set_xlim(x_start, x_end)
ax.set_xscale("log", basex=2)
ax.set_yscale("log", basey=2)
ax.set_facecolor('#E6E6E6')
ax.grid(False)
ax.grid(which='major', axis='y', color='w', linestyle='-')


# plt.legend()

#alternative directly in plot
ax.annotate('Peak Performance \u03A0 (vector)', xy=(32, system_info_dict['peak_perf_vec']), xytext=(32, system_info_dict['peak_perf_vec']+2), fontsize=10, color='forestgreen')
ax.annotate('Peak Performance \u03A0 (scalar)', xy=(4, system_info_dict['peak_perf']), xytext=(4, system_info_dict['peak_perf']+0.2), fontsize=10, color='mediumseagreen')

ax.annotate('Read bandwidth theoretical \u03B2', xy=(0.5, system_info_dict['mem_bandwidth_theoretical']*0.5), xytext=(0.3, system_info_dict['mem_bandwidth_theoretical']*0.4), fontsize=10, color='maroon', rotation=58)
ax.annotate('Read bandwidth measured \u03B2', xy=(0.5, system_info_dict['mem_bandwidth_measured']*0.5), xytext=(0.3, system_info_dict['mem_bandwidth_measured']*0.4), fontsize=10, color='indianred', rotation=58)


plt.subplots_adjust(top=0.78)

plt.tight_layout()

fig = plt.gcf()
plt.show()
# fig.savefig("test/benchmark_plot.png", bbox_inches="tight")
fig.savefig("test/benchmark_plot.pdf", bbox_inches="tight")
print("----DONE----")