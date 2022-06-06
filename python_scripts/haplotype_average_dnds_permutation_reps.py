#!/usr/bin/python3

import argparse
import os
import sys
from functions.io_utils import read_fasta
from functions.dn_ds import pairwise_dnds
import itertools
import pandas as pd
import numpy as np

def main():

    parser = argparse.ArgumentParser(

    description="Same as other script, but will run 1000 replicates based on permuting the haplotype numbers across the abundance table. Otherwise: calculates dN/dS averaged across subset of haplotypes in each sample with at least two haplotypes). "
                "Will also compute dN/dS for haplotypes within and not within each given sample. ",

    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input_haplotypes', metavar='FASTA', type=str,
                        help="Path to codon-aligned haplotype FASTA.", required=True)

    parser.add_argument('-a', '--abun', metavar='TABLE', type=str,
                        help="Path to haplotype abundance table (across samples).",
                        required=True)

    parser.add_argument('-s', '--samples', metavar='SAMPLES', type=str,
                        help="Path to file of samples used with StrainFinder.",
                        required=True)

    parser.add_argument('-g', '--gene', metavar='GENE', type=str,
                        help='Gene name (used for matching haplotype number to haplotype name in input FASTA).', required=True)

    parser.add_argument('-m', '--min_cutoff', metavar='FLOAT', type=float,
                        help='Min abundance cut-off for calling haplotypes as present',
                        required=True)

    parser.add_argument("-o", "--out_prefix", metavar="OUTPUT", type=str,
                        help="Path to output prefix.", required=True)
    
    args = parser.parse_args()


    samples = []
    with open(args.samples, 'r') as sample_fh:
        for sample_line in sample_fh:
            samples.append(sample_line.rstrip())

    seqs = read_fasta(args.input_haplotypes)

    if len(seqs.keys()) == 1:
        sys.exit("Stopping as there is only one gene haplotype present.")

    pairwise_combos = list(itertools.combinations(seqs, 2))

    # Get table of all pairwise dn/ds values.
    haplotype_dnds = pd.DataFrame(columns = ['dn', 'ds', 'dnds'],
                                  index = pd.MultiIndex.from_tuples(pairwise_combos,
                                                                    names=('seq1', 'seq2')))

    haplotype_indices_level0 = list(haplotype_dnds.index.get_level_values(0))

    haplotype_indices_level1 = list(haplotype_dnds.index.get_level_values(1))

    for combo in pairwise_combos:

        haplotype_dnds.loc[combo, :] = list(pairwise_dnds(seqs[combo[0]],
                                                          seqs[combo[1]]))

    haplotype_abun = pd.read_csv(args.abun, sep = "\t")

    haplotype_abun[haplotype_abun < args.min_cutoff] = 0

    haplotype_abun.index = samples

    all_mean_per_sample_dnds = pd.DataFrame(columns = ['within_dnds', 'between_dnds'],
                                            index = list(range(1, 1001)))

    for rep_i in range(1, 1001):

        per_sample_dnds = pd.DataFrame(columns = ['within_dnds', 'between_dnds'], index = samples)

        rep_ran_columns = list(haplotype_abun.columns)
        np.random.shuffle(rep_ran_columns)

        ran_haplotype_abun = haplotype_abun.loc[:, rep_ran_columns]

        for s in samples:
            
            sample_haplotype_abun = ran_haplotype_abun.loc[s, :]

            haplotypes_present = list(np.nonzero(list(sample_haplotype_abun))[0])

            if len(haplotypes_present) > 1:

                haplotypes_present_ids = []

                for haplotype_num in haplotypes_present:

                    haplotype_num = haplotype_num + 1

                    haplotypes_present_ids.append("haplotype" + str(haplotype_num) + "_" + args.gene)

                indices_with_present_haplotypes = []
                indices_with_only_one_present_haplotype = []

                for i in range(len(haplotype_indices_level0)):

                    if haplotype_indices_level0[i] in haplotypes_present_ids and haplotype_indices_level1[i] in haplotypes_present_ids:
                        indices_with_present_haplotypes.append(i)
                    elif haplotype_indices_level0[i] in haplotypes_present_ids and haplotype_indices_level1[i] not in haplotypes_present_ids:
                        indices_with_only_one_present_haplotype.append(i)
                    elif haplotype_indices_level0[i] not in haplotypes_present_ids and haplotype_indices_level1[i] in haplotypes_present_ids:
                        indices_with_only_one_present_haplotype.append(i)

                if len(indices_with_present_haplotypes) > 1 and haplotype_dnds.iloc[indices_with_present_haplotypes, 2].count() > 0:
                    per_sample_dnds.loc[s, 'within_dnds'] = np.nanmean(np.array(haplotype_dnds.iloc[indices_with_present_haplotypes, 2]), dtype='float32')

                if len(indices_with_only_one_present_haplotype) > 1 and haplotype_dnds.iloc[indices_with_only_one_present_haplotype, 2].count() > 0:
                    per_sample_dnds.loc[s, 'between_dnds'] = np.nanmean(np.array(haplotype_dnds.iloc[indices_with_only_one_present_haplotype, 2]), dtype='float32')

        if per_sample_dnds.within_dnds.count() > 0:
            all_mean_per_sample_dnds.loc[rep_i, 'within_dnds'] = np.nanmean(np.array(per_sample_dnds.within_dnds), dtype='float32')

        if per_sample_dnds.between_dnds.count() > 0:
            all_mean_per_sample_dnds.loc[rep_i, 'between_dnds'] = np.nanmean(np.array(per_sample_dnds.between_dnds), dtype='float32')

    # Write output.
    all_mean_per_sample_dnds.to_csv(args.out_prefix + "_rep_mean_dnds.tsv", sep = "\t",
                                    header = True, index_label = 'replicate', na_rep='NA')

if __name__ == '__main__':
    main()
