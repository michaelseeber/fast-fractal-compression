import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import time
import re
import json
import seaborn as sns

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

labels = ("Manual vectorization + entropy: "+ r"$\bf{0.59s}$", "Manual vectorization: "+ r"$\bf{4.7s}$","Fast unoptimized: "+ r"$\bf{40.2s}$", "Base: "+ r"$\bf{229.4s}$")
y_pos = np.arange(len(labels))
runtime = [0.59,4.7,40.2,229.4]

sns.set()
sns.set_palette(sns.hls_palette(8, l=.3, s=.8), color_codes=True)

fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(1, 1, 1)

plt.title(SYSTEM_INFO, loc="left", pad=10, fontweight="bold")
# ax.set_xlabel('Image size [Pixel]')
# ax.set_ylabel('Performance [Flops/Cycle]', rotation=0)

plt.barh(y_pos, runtime, align='center', alpha=1, color=['firebrick', 'navy', 'olivedrab', 'teal'])
plt.yticks([])
plt.xlabel('Compression time [seconds]')


# # manual annotation of names 
plt.annotate(labels[3], xy=(95 , 2.9), color='white', fontsize="14")
plt.annotate(labels[2], xy=(43 , 1.9), color='black',fontsize="14")
plt.annotate(labels[1], xy=(7 , 0.9), color='black', fontsize="14")
plt.annotate(labels[0], xy=(3 , -0.1), color='black', fontsize="14")
# ax.annotate(labels[3], xy=(1000 , 10.5), color='olivedrab')
# ax.annotate(labels[4], xy=(1400, 13), color='darkgoldenrod')


# ax.set_xlim(100,2100)
# ax.set_xticks([256, 512, 1024, 2048])
# ax.set_ylim(0,20)
ax.set_facecolor('#E6E6E6')
# ax.grid(False)
# ax.grid(which='major', axis='y', color='w', linestyle='-')
# plt.legend()

plt.tight_layout()
fig = plt.gcf()
plt.show()
fig.savefig("test/runtime_plot_report.pdf",  bbox_inches="tight")
# fig.savefig("test/performance_plot.png", bbox_inches="tight")
print("----DONE----")
