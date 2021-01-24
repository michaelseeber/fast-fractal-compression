import numpy as np
import matplotlib.pyplot as plt
import subprocess
import seaborn as sns
import time
import re
import json

from subprocess import Popen, PIPE


def run_fic(version, input_img):

    #p = Popen([version["cmd"], "images/" + input_img], cwd=working_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    p = Popen([version["cmd"], "images/report/" + input_img], cwd=working_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    print("Running implementation " + version["name"] +" with input " + input_img +" ...")
    fic_output, fic_errors = p.communicate()

    if fic_errors != '':
        print(fic_errors)
    else:
        return fic_output



#Constants and arguments
try:
    with open('test/system_info.json') as f:
        system_info_dict = json.load(f)
except Exception as e:
    print(e)
    print('')
    print('Please create your system info file. You can use test/system_info_michi.json as an example.')
    exit(1)
SYSTEM_INFO = "" + system_info_dict['cpu_name'] + ' @ ' + system_info_dict['cpu_base_freq'] \
    + '\nCompiler: ' + system_info_dict['compiler'] + ', OS: ' + system_info_dict['os']
#fic_reference = { "name" : "fic_reference", "cmd" : "./release/fic_compress_flip"}
#fic_fast= { "name": "fic_fast", "cmd" : "./release/fic_compress_flip_fast"}
fic_fast_reference = { "name": "fic_fast_reference", "cmd" : "./release_icc/fic_compress_flip_fast_reference"}
fic_scalar_replacement = { "name": "fic_scalar_replacement", "cmd" : "./release/fic_compress_flip_fast_scalar_replacement"}
fic_inline = { "name": "fic_inline", "cmd" : "./release/fic_compress_flip_fast_inline"}
fic_ilp = { "name": "fic_ilp", "cmd" : "./release/fic_compress_flip_fast_ilp"}
fic_avx = { "name": "fic_avx", "cmd" : "./release/fic_compress_flip_fast_avx"}
compile_cmd = "make"
working_dir_compile = 'fic/release'
working_dir = 'fic/'

#input_images = [{ "name": "256_aerial.jpg", "size" : 256 },
#                { "name": "512_bridge.jpg", "size" :512 },
#                { "name": "1024_airport.jpg", "size" : 1024 },
#                { "name": "2048_road.jpg", "size" : 2048 }]
input_images = [{ "name": "128_ape.bmp", "size" : 128 },
                { "name": "256_ape.bmp", "size" : 256 },
                { "name": "384_ape.bmp", "size" : 384 },
                { "name": "512_ape.bmp", "size" : 512 },
                { "name": "640_ape.bmp", "size" : 640 },
                { "name": "768_ape.bmp", "size" : 768 },
                { "name": "896_ape.bmp", "size" : 896 },
                { "name": "1024_ape.bmp", "size" : 1024 },
                { "name": "1152_ape.bmp", "size" : 1152 },
                { "name": "1280_ape.bmp", "size" : 1280 },
                { "name": "1408_ape.bmp", "size" : 1408 },
                { "name": "1536_ape.bmp", "size" : 1536 },
                { "name": "1664_ape.bmp", "size" : 1664 },
                { "name": "1792_ape.bmp", "size" : 1792 },
                { "name": "1920_ape.bmp", "size" : 1920 },
                { "name": "2048_ape.bmp", "size" : 2048 }]

#load from log or recompile
answer = input("Yes to compile and run. No to plot form logs. (y/n) ")
if answer == "y":
    #compile code
    print("Compiling...")
    p = Popen(compile_cmd, cwd=working_dir_compile, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    compile_output, compile_errors = p.communicate()

    if compile_errors != '':
        print(compile_errors)
    else:
        print(compile_output)

    #run program
    print("Running...")
    fic_result = ""
    for img in input_images:
        fic_result += "Image: " + img["name"] + "/n"
        #fic_result += run_fic(fic_reference, img["name"])
        #fic_result += run_fic(fic_fast, img["name"])
        fic_result += run_fic(fic_fast_reference, img["name"])
        fic_result += run_fic(fic_scalar_replacement, img["name"])
        fic_result += run_fic(fic_inline, img["name"])
        fic_result += run_fic(fic_ilp, img["name"])
        fic_result += run_fic(fic_avx, img["name"])

    #store for logging
    path_to_output_file = 'test/results_performance/result_' + time.strftime("%d.%m_%H:%M") +'.txt'
    with open(path_to_output_file, 'w+') as output_file:
        output_file.write(fic_result)
elif answer == "n":

    report = input("You want to create the report plot?(y/n) ")

    if report == "y":
        path_to_log_file = 'test/results_performance/result_merged.txt'

    elif report == "n":       
        log_date = input("Enter time of log. Example: 26.04_15:40  ")
        path_to_log_file = 'test/results_performance/result_' + log_date +'.txt'

    else:
        print("No valid input")
        exit(1)

    try:
        with open(path_to_log_file, 'r') as log:
            fic_result = log.read()
    except IOError:
        print ("Could not read file: " + path_to_log_file)
        exit(1)
else:
    print("No valid input")
    exit(1)



#plotting code
#labels = ["base_flip","fast_flip"]
labels = ["Fast unoptimized", "Scalar replacement","Inlining", "ILP", "Manual vectorization"]
parse = re.findall("Running:\s*(\S*)\s(\d*.?\d*)", fic_result)

sns.set()
sns.set_palette(sns.hls_palette(8, l=.3, s=.8), color_codes=True)

fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(1, 1, 1)

ax.set_title(SYSTEM_INFO, loc="left", pad=30, fontweight="bold")
ax.set_xlabel('Image size [Pixel]')
ax.set_ylabel('Performance [Flops/Cycle]', rotation=0)
ax.yaxis.set_label_coords(0.15, 1.02)

#ax.plot([el["size"] for el in input_images], [float(el[1]) for el in parse[::2]], 'o-', label=labels[0], color='firebrick')
#ax.plot([el["size"] for el in input_images], [float(el[1]) for el in parse[1::2]], 'o-', label=labels[1], color='navy')
ax.plot([el["size"] for el in input_images], [float(el[1]) for el in parse[::5]], 'o-', label=labels[0], color='firebrick')
ax.plot([el["size"] for el in input_images], [float(el[1]) for el in parse[1::5]], 'o-', label=labels[1], color='navy')
ax.plot([el["size"] for el in input_images], [float(el[1]) for el in parse[2::5]], 'o-', label=labels[2], color='teal')
ax.plot([el["size"] for el in input_images], [float(el[1]) for el in parse[3::5]], 'o-', label=labels[3], color='olivedrab')
ax.plot([el["size"] for el in input_images], [float(el[1]) for el in parse[4::5]], 'o-', label=labels[4], color='darkgoldenrod')

# manual annotation of names 
ax.annotate(labels[0], xy=(1100 , 2.25), color='firebrick')
ax.annotate(labels[1], xy=(500 , 3.5), color='navy')
ax.annotate(labels[2], xy=(1200 , 7), color='teal')
ax.annotate(labels[3], xy=(1000 , 10.5), color='olivedrab')
ax.annotate(labels[4], xy=(1300, 12), color='darkgoldenrod')
#for icc
# ax.annotate(labels[4], xy=(1400, 11.5), color='darkgoldenrod')

ax.set_xlim(100,2100)
ax.set_xticks([256, 512, 1024, 2048])
ax.set_ylim(0,20)
ax.set_facecolor('#E6E6E6')
ax.grid(False)
ax.grid(which='major', axis='y', color='w', linestyle='-')
# plt.legend()


plt.subplots_adjust(top=0.78)

plt.tight_layout()
fig = plt.gcf()
plt.show()
# fig.savefig("test/performance_plot_report.pdf", bbox_inches="tight")
# fig.savefig("test/performance_plot_icc_report.pdf", bbox_inches="tight")
# fig.savefig("test/performance_plot.png", bbox_inches="tight")
print("----DONE----")
