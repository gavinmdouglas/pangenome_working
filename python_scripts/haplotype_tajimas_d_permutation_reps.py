#!/usr/bin/python3

import argparse
import os
import sys
from functions.io_utils import read_fasta, write_fasta
from functions.pop_gen import tajimas_d_and_diversity_metrics
import pandas as pd
import numpy as np


def main():

    parser = argparse.ArgumentParser(

    description="Same as other script, but based on permuting allele columns for 1000 replicates and returning output. "
                "Calc Tajima's D and nucleotide diversity metrics for all samples with at least two haplotypes. "
                "Takes mean abundance of genes per sample into account to get \"counts\" of "
                "abundance based on the StrainFinder-outputted relative abundance.",

    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input_haplotypes', metavar='FASTA', type=str,
                        help='Path to codon-aligned haplotype FASTA.', required=True)

    parser.add_argument('-a', '--abun', metavar='TABLE', type=str,
                        help='Path to haplotype abundance table (across samples).',
                        required=True)

    parser.add_argument('-d', '--depth', metavar='TABLE', type=str,
                        help='Path to gene mean depth table (genes as rows, samples as columns)',
                        required = True)

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

    # Now get mean gene depth (based on variant sites) per sample.
    all_mean_depth = pd.read_csv(filepath_or_buffer = args.depth,
                                 sep = '\t', index_col = 0)

    gene_mean_depth = all_mean_depth.loc[args.gene, samples]

    haplotype_abun = pd.read_csv(args.abun, sep = "\t")

    haplotype_abun[haplotype_abun < args.min_cutoff] = 0

    haplotype_abun.index = samples

    all_mean_per_sample_tajimas_D = pd.DataFrame(columns = ['mean_S', 'mean_Wattersons', 'mean_pi', 'mean_D'],
                                                 index = list(range(1, 1001)))

    for rep_i in range(1, 1001):

        rep_ran_columns = list(haplotype_abun.columns)
        np.random.shuffle(rep_ran_columns)

        ran_haplotype_abun = haplotype_abun.loc[:, rep_ran_columns]

        # Get table of all diversity metrics.
        haplotype_div_metrics = pd.DataFrame(columns = ['n', 'S', 'Wattersons', 'pi','D'],
                                             index = samples)

        for s in samples:
            
            sample_haplotype_abun = ran_haplotype_abun.loc[s, :]

            haplotypes_present = list(np.nonzero(list(sample_haplotype_abun))[0])

            if len(haplotypes_present) > 1:

                haplotypes_present_seqs = []

                for haplotype_num in haplotypes_present:

                    haplotype_relabun = sample_haplotype_abun[haplotype_num]

                    haplotype_depth = round(haplotype_relabun * gene_mean_depth[s])

                    # Set minimum depth to be 1.
                    if haplotype_depth < 1:
                        haplotype_depth = 1

                    haplotype_num = haplotype_num + 1

                    # Add in haplotype sequence the number of times equal to the estimated depth.
                    for i in range(haplotype_depth):
                        haplotypes_present_seqs.append(seqs["haplotype" + str(haplotype_num) + "_" + args.gene])

                haplotype_div_metrics.loc[s, :] = list(tajimas_d_and_diversity_metrics(haplotypes_present_seqs))

            else:
                haplotype_div_metrics.loc[s, :] = [len(haplotypes_present), float("NaN"), float("NaN"), float("NaN"), float("NaN")]

        if haplotype_div_metrics.S.count() > 0:
            all_mean_per_sample_tajimas_D.loc[rep_i, 'mean_S'] = np.nanmean(np.array(haplotype_div_metrics.S), dtype='float32')

        if haplotype_div_metrics.Wattersons.count() > 0:
            all_mean_per_sample_tajimas_D.loc[rep_i, 'mean_Wattersons'] = np.nanmean(np.array(haplotype_div_metrics.Wattersons), dtype='float32')

        if haplotype_div_metrics.pi.count() > 0:
            all_mean_per_sample_tajimas_D.loc[rep_i, 'mean_pi'] = np.nanmean(np.array(haplotype_div_metrics.pi), dtype='float32')

        if haplotype_div_metrics.D.count() > 0:
            all_mean_per_sample_tajimas_D.loc[rep_i, 'mean_D'] = np.nanmean(np.array(haplotype_div_metrics.D), dtype='float32')

    # Write output.
    all_mean_per_sample_tajimas_D.to_csv(args.out_prefix + "_rep_mean_results.tsv",
                                         sep = "\t", header = True, na_rep='NA',
                                         index_label = 'category')

if __name__ == '__main__':
    main()
