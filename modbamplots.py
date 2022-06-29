#! /usr/bin/env python

from modbampy import ModBam
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--bam', help='the modified bam file that has been generated', required=True)
parser.add_argument('--chr', help='the chromosome', required=True)
parser.add_argument('--start', help='the start coordinates', required=True)
parser.add_argument('--end', help='the end coordinates', required=True)

args = parser.parse_args()

with ModBam(args.bam) as bam:
    for read in bam.reads(args.chrom, args.start, args.end):
        for pos_mod in read.mod_sites:
            print(*pos_mod)
