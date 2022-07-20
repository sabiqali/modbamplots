#! /usr/bin/env python

from modbampy import ModBam
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

parser = argparse.ArgumentParser()
parser.add_argument('--bam', help='the modified bam file that has been generated', required=False)
parser.add_argument('--chr', help='the chromosome', required=False)
parser.add_argument('--start', help='the start coordinates', required=False)
parser.add_argument('--end', help='the end coordinates', required=False)

args = parser.parse_args()

with ModBam(args.bam) as bam:
    bam_meth = {} #{pos : [min,avg,max,read_support]}
    for read in bam.reads(args.chr, int(args.start), int(args.end)):
        for pos_mod in read.mod_sites:
            #print(*pos_mod)
            pos_count = pos_mod.rpos
            inst_meth = (pos_mod.qual)/255
            #print(pos_count)
            #print(inst_meth)
            if pos_count > 0:
                if pos_count in bam_meth:
                    #update that by averaging and checking min max
                    pos_vals = bam_meth[pos_count]
                    new_min = inst_meth if inst_meth < pos_vals[0] else pos_vals[0]
                    new_max = inst_meth if inst_meth > pos_vals[2] else pos_vals[2]
                    new_avg = ((pos_vals[1] * pos_vals[3]) + inst_meth)/(pos_vals[3] + 1)
                    new_read_support = pos_vals[3] + 1
                    bam_meth.update({pos_count : [new_min, new_avg, new_max, new_read_support]})
                else:
                    bam_meth[pos_count] = [inst_meth, inst_meth, inst_meth, 1]

#print(bam_meth)
#print(sorted(bam_meth.items()))
#sorted_bam_meth = sorted(bam_meth.items())
bam_meth_keys = list(bam_meth.keys())
bam_meth_keys.sort()
time = list()
minimum = list()
maximum = list()
average = list()
support = list()
for key in bam_meth_keys:
    #print(key)
    #print(bam_meth[key][1])
    if key >= int(args.start) and key <= int(args.end):
        time.append(key)
        minimum.append(bam_meth[key][0])
        average.append(bam_meth[key][1])
        maximum.append(bam_meth[key][2])
        support.append(bam_meth[key][3])
#time = np.arange(int(args.end) - int(args.start))
#time2 = np.array([1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011])
#min = np.array([1, 3, 1, 4, 5, 2, 4, 3, 2, 5, 1, 5])
#min = np.random.rand(int(args.end) - int(args.start))
#max = np.random.rand(int(args.end) - int(args.start))
#avg = np.random.rand(int(args.end) - int(args.start))

# Initialize figure and axis
fig, (ax1,ax2) = plt.subplots(2,1,sharex=True,figsize=(8, 8))

# Plot lines
ax1.plot(time, minimum, color="red")
ax1.plot(time, maximum, color="red")
ax1.plot(time, average, color="black")
ax1.fill_between(time, minimum, maximum, color='red',alpha=0.3,label="Variance in Methylation")

#ax1.ylabel("methylation")

ax2.plot(time, support, color="green")

#ax2.ylabel("read support")

ax1.set_title('Methylation')
ax2.set_title('Read Support')

#plt.fill_between(time, minimum, maximum, color='red',alpha=0.3,label="Variance in Methylation")
# Fill area when income > expenses with green
#ax.fill_between(
#    time, minimum, maximum, where=(maximum > minimum), 
#    interpolate=True, color="red", alpha=0.25, 
#    label="Variance in Methylation"
#)
# Creating legend with color box
red_patch = mpatches.Patch(color='red', label='Min & Max Methylation')
black_patch = mpatches.Patch(color='black', label='Average Methylation')
green_patch = mpatches.Patch(color='green', label='Read Support')
plt.legend(handles=[red_patch, black_patch, green_patch])


plt.xlabel("Region")

fig.savefig("test.png")
