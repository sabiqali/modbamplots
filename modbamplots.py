#! /usr/bin/env python

from modbampy import ModBam
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np

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
            print(*pos_mod)
            pos_count = pos_mod.rpos
            inst_meth = pos_mod.qual
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
        break

time = np.arange(int(args.end) - int(args.start))
#time2 = np.array([1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011])
#min = np.array([1, 3, 1, 4, 5, 2, 4, 3, 2, 5, 1, 5])
min = np.random.rand(int(args.end) - int(args.start))
max = np.random.rand(int(args.end) - int(args.start))
avg = np.random.rand(int(args.end) - int(args.start))

# Initialize figure and axis
fig, ax = plt.subplots(figsize=(8, 8))

# Plot lines
ax.plot(time, min, color="red")
ax.plot(time, max, color="red")
ax.plot(time, avg, color="black")

# Fill area when income > expenses with green
ax.fill_between(
    time, min, max, where=(max > min), 
    interpolate=True, color="red", alpha=0.25, 
    label="Variance in Methylation"
)

plt.xlabel("Region")
plt.ylabel("Methylation")

ax.legend()

fig.savefig("test.png")
