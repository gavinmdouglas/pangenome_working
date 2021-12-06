#!/usr/bin/python3

import argparse
import os
import sys
from functions.io_utils import read_fasta
from functions.codon import start_codon_present
import pandas as pd
import numpy as np


def main():

    parser = argparse.ArgumentParser(

    description='Call start / stop codon knockout mutations and premature stop '
                'codons in gene FASTAs. Assumes that each input sequence '
                '*should* start with start codon and end in stop codon. Note '
                'that start and stop codons will only be called as missing if '
                'there is not another start/stop within the neighbouring three '
                'codons.',

    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input', metavar='FASTA', type=str,
                        help='DNA sequence FASTA of input genes.',
                        required=True)

    parser.add_argument('-o', '--output', metavar='TABLE', type=str,
                        help='Output table to create.', required=True)

    args = parser.parse_args()

    seqs = read_fasta(args.input)

    seq_info = pd.DataFrame(columns = ['start_codon_missing', 'stop_codon_missing', 'premature_stop_codon_number', 'expected_stop_codon_number', 'percent_truncated'],
                            index = seqs.keys())
       

    for seq_id, seq in seqs.items():

        seq = seq.replace('-', '')

        seq = seq.upper()

        (start_codon_present, start_codon_position) = start_codon_present(seq)

        start_codon_missing = not start_codon_present




        print(start_codon_missing)


   #haplotype_dnds.to_csv(args.out_prefix + "_pairwise_haplotype_dnds.tsv",
    #                      sep = "\t", header = True, na_rep='NA')

if __name__ == '__main__':
    main()
