#!/usr/bin/python3

import argparse
import os
import sys
from functions.io_utils import read_fasta
from functions.pop_gen import tajimas_d_and_diversity_metrics
import pandas as pd
import numpy as np


def main():

    parser = argparse.ArgumentParser(

    description="Calc Tajima's D and nucleotide diversity metrics (and also for the subset of haplotypes in each sample with at least three haplotypes). Also calculate based on all gene sequences from reference genomes.",

    epilog='''Usage example:

    python haplotype_tajimas_d.py 

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input_haplotypes', metavar='FASTA', type=str,
                        help='Path to codon-aligned haplotype FASTA.', required=True)

    parser.add_argument('-a', '--abun', metavar='TABLE', type=str,
                        help='Path to haplotype abundance table (across samples).',
                        required=True)

    parser.add_argument('-r', '--ref', metavar='FASTA', type=str,
                        help='Path to codon-aligned FASTA for genes from reference genome.',
                        required=True)  

    parser.add_argument('-s', '--samples', metavar='SAMPLES', type=str,
                        help='Path to file of samples used with StrainFinder.',
                        required=True)

    parser.add_argument('-g', '--gene', metavar='GENE', type=str,
                        help='Gene name (used for matching haplotype number to haplotype name in input FASTA).', required=True)

    parser.add_argument('-m', '--min_cutoff', metavar='FLOAT', type=float,
                        help='Min abundance cut-off for calling haplotypes as present',
                        required=True)

    parser.add_argument('-o', '--out_prefix', metavar='OUTPUT', type=str,
                        help='Path to output prefix.', required=True)
    
    args = parser.parse_args()

    samples = []
    with open(args.samples, 'r') as sample_fh:
        for sample_line in sample_fh:
            samples.append(sample_line.rstrip())

    seqs = read_fasta(args.input_haplotypes)
    ref_seqs = read_fasta(args.ref)

    #Get table of all diversity metrics.
    haplotype_div_metrics = pd.DataFrame(columns = ['n', 'S', 'Wattersons', 'pi','D'],
                                        index = samples + ['all_haplotypes',
                                                           'reference_seqs'])

    haplotype_div_metrics.loc['all_haplotypes', :] = list(tajimas_d_and_diversity_metrics(seqs.values()))

    sys.exit()
    haplotype_div_metrics.loc['reference_seqs', :] = list(tajimas_d_and_diversity_metrics(ref_seqs.values()))


    haplotype_abun = pd.read_csv(args.abun, sep = "\t")

    haplotype_abun[haplotype_abun < args.min_cutoff] = 0

    haplotype_abun.index = samples

    for s in samples:
        
        sample_haplotype_abun = haplotype_abun.loc[s, :]

        haplotypes_present = list(np.nonzero(list(sample_haplotype_abun))[0])

        if len(haplotypes_present) > 2:

            haplotypes_present_seqs = []

            for haplotype_num in haplotypes_present:

                haplotype_num = haplotype_num + 1

                haplotypes_present_seqs.append(seqs["haplotype" + str(haplotype_num) + "_" + args.gene])

            haplotype_div_metrics.loc[s, :] = list(tajimas_d_and_diversity_metrics(haplotypes_present_seqs))

        else:

            haplotype_div_metrics.loc[s, :] = [len(haplotypes_present), float("NaN"), float("NaN"), float("NaN"), float("NaN")]


    # Write output.
    haplotype_div_metrics.to_csv(args.out_prefix + "_tajimas_d_and_metrics.tsv",
                                 sep = "\t", header = True, na_rep='NA',
                                 index_label = 'category')

if __name__ == '__main__':
    main()
