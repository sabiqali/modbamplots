#! /usr/bin/env python

from modbampy import ModBam
import argparse

parser = argparse.ArgumentParser(description='Find out the methylation stats of a Protein Coding Region from BAM files with Modified Base tags, and then plot methylation by distance in the region')
parser.add_argument('--bam', help='the modified bam file that has been generated', required=True)
parser.add_argument('--region', help='a tsv containing the protein coding gene details', required=True)

args = parser.parse_args()

region_fh = open(args.region)

header = region_fh.readline()

print("\t".join(["chromosome","gene","absolute_position","region","relative_position","relative_scaled","methylation_probability","read_support"]))

for line in region_fh:
    gene_name,chromosome,upstream_start,tss,tes,downstream_end = line.rstrip().split('\t')
    with ModBam(args.bam) as bam:
        bam_meth = {} #{pos : [min,avg,max,read_support]}
        for read in bam.reads(chromosome, int(upstream_start), int(downstream_end)):
            for pos_mod in read.mod_sites:
                pos_count = pos_mod.rpos
                inst_meth = (pos_mod.qual)/255
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

    for position in bam_meth:
        if int(position) >= int(upstream_start) and int(position) < int(tss) :
            region_of_mod = "upstream"
            rel_region_of_mod = int(position) - int(upstream_start)
            rel_region_of_mod_scaled = (int(position) - int(upstream_start))/(int(tss) - int(upstream_start))
            meth_prob = bam_meth[position][1]
            read_support = bam_meth[position][3]
        if int(position) >= int(tss) and int(position) <= int(tes) :
            region_of_mod = "transcript"
            rel_region_of_mod = int(position) - int(tss)
            rel_region_of_mod_scaled = (int(position) - int(tss))/(int(tes) - int(tss))
            meth_prob = bam_meth[position][1]
            read_support = bam_meth[position][3]
        if int(position) > int(tes) and int(position) <= int(downstream_end) :
            region_of_mod = "downstream"
            rel_region_of_mod = int(position) - int(tes)
            rel_region_of_mod_scaled = (int(position) - int(tes))/(int(downstream_end) - int(tes))
            meth_prob = bam_meth[position][1]
            read_support = bam_meth[position][3]
        
        print("\t".join([chromosome,gene_name,str(position),region_of_mod,rel_region_of_mod,rel_region_of_mod_scaled,meth_prob,read_support]))

    break

#print(bam_meth)