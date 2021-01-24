import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import time
import re
import json

from subprocess import Popen, PIPE


def run_fic(version, input_img):

    p = Popen([version["cmd"], "images/" + input_img], cwd=working_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
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
fic_reference = { "name" : "fic_reference", "cmd" : "./release/fic_compress"}
fic_fast= { "name": "fic_fast", "cmd" : "./release/fic_compress_fast"}
compile_cmd = "make"
working_dir_compile = 'fic/release'
working_dir = 'fic/'

input_images = [{ "name": "256_aerial.jpg", "size" : 256 },
                { "name": "512_bridge.jpg", "size" :512 },
                { "name": "1024_airport.jpg", "size" : 1024 },
                { "name": "2048_road.jpg", "size" : 2048 }]

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
        fic_result += run_fic(fic_reference, img["name"])
        fic_result += run_fic(fic_fast, img["name"])

    #store for logging
    path_to_output_file = 'test/results_performance/result_' + time.strftime("%d.%m_%H:%M") +'.txt'
    with open(path_to_output_file, 'w+') as output_file:
        output_file.write(fic_result)
elif answer == "n":
    log_date = input("Enter time of log. Example: 26.04_15:40  ")
    path_to_log_file = 'test/results_performance/result_' + log_date +'.txt'
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
labels = ["base","fast"]
parse = re.findall("Running:\s*(\S*)\s(\d*.?\d*)", fic_result)


sns.set()
sns.set_palette(sns.hls_palette(8, l=.3, s=.8), color_codes=True)

fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(1, 1, 1)

ax.set_title(SYSTEM_INFO, loc="left", pad=30, fontweight="bold")
ax.set_xlabel('Image size [Pixel]')
ax.set_ylabel('Performance [Flops/Cycle]', rotation=0)
ax.yaxis.set_label_coords(0.15, 1.02)


ax.plot([el["size"] for el in input_images], [float(el[1]) for el in parse[::2]], 'o-', label=labels[0], color='firebrick')
ax.plot([el["size"] for el in input_images], [float(el[1]) for el in parse[1::2]], 'o-', label=labels[1], color='navy')


ax.set_xlim(200,2200)
ax.set_xticks([256, 512, 1024, 2048])
ax.set_ylim(0,16)
ax.set_facecolor('#E6E6E6')
ax.grid(False)
ax.grid(which='major', axis='y', color='w', linestyle='-')
plt.legend()
plt.subplots_adjust(top=0.78)

plt.tight_layout()
fig = plt.gcf()
plt.show()
fig.savefig("test/performance_plot.png", bbox_inches="tight")
fig.savefig("test/performance_plot.pdf", bbox_inches="tight")
print("----DONE----")
