#!/usr/bin/python3

import argparse
import os
import sys
import pysam

bam_path = "/home/gdouglas/projects/honey_bee/Ellegaard/pangenome_competitive_mapping/bowtie2_mapping/panaroo/mapping/SRR7287219.Gilliamella.merged.nonparalog.bam"

bam = pysam.AlignmentFile(bam_path, 'rb')

contig_name = "Gilli_group_1710"
gene_start_coor = 100
gene_end_coor = 359

gene_length = gene_end_coor - gene_start_coor
num_covered_sites = 0

for pileupcolumn in bam.pileup(contig_name, gene_start_coor, gene_end_coor, stepper = 'nofilter', truncate = True):

    # Skip any positions not in the gene range.
    if pileupcolumn.pos < gene_start_coor or pileupcolumn.pos > gene_end_coor:
        continue
 
    pos_coverage = len(pileupcolumn.pileups)

    if pos_coverage > 0:
        num_covered_sites += 1

    print(str(pileupcolumn.pos) + "\t" + str(pos_coverage))


print("\t\t\t")
print("breadth: " + str(num_covered_sites / gene_length))
