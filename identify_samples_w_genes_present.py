#!/usr/bin/python3

import argparse
import os
import pandas as pd

def main():

    parser = argparse.ArgumentParser(

            description="Parse gene ids and PANDORA matrix and figure out in which samples a minimum prop of the genes were called as present.",

    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-p", "--presence_matrix", metavar="MATRIX", type=str,
                        help="Path to gene presence absence matrix (samples are columns)", required=True)

    parser.add_argument("-g", "--gene_ids", metavar="GENE_IDS", type=str,
                        help="Path to file with gene ids.",
                        required=True)

    parser.add_argument("-c", "--cutoff", metavar="CUTOFF", type=float,
                        help="Prop of genes that must be present for a sample to be retained.",
                        required=False, default = 0.95)

    parser.add_argument("-o", "--output", metavar="OUTPUT", type=str,
                        help="Path to output file with samples that ",
                        required=True)

    args = parser.parse_args()

    genes_to_keep = []
    with open(args.gene_ids, 'r') as GENE_IDS:
        for line in GENE_IDS:
            genes_to_keep.append(line.rstrip())

    presence_tab = pd.read_csv(args.presence_matrix, sep = "\t", index_col = 0)

    num_genes_present = presence_tab.loc[genes_to_keep, :].sum(axis = 0)

    min_num_genes = len(genes_to_keep) * args.cutoff

    samples_to_keep = list(num_genes_present.loc[num_genes_present >= min_num_genes].index)

    with open(args.output, "w") as OUTPUT:
        for s in samples_to_keep:
            print(s, file = OUTPUT)


if __name__ == '__main__':

    main()
