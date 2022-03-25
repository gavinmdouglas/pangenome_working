#!/usr/bin/python3

import argparse
import os
import sys
from functions.io_utils import read_fasta
from functions.codon import start_codon_present, stop_codon_premature_present
import pandas as pd
import numpy as np
import math

def main():

    parser = argparse.ArgumentParser(

    description='Call canonical start / stop codon knockout mutations and premature stop '
                'codons in gene FASTAs. Assumes that each input sequence '
                '*should* start with start codon and end in stop codon.',

    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input', metavar='FASTA', type=str,
                        help='DNA sequence FASTA of input genes.',
                        required=True)

    parser.add_argument('-n', '--name', metavar='NAME', type=str,
                        help='Name to include in output table: usually the unique gene or ortholog id.',
                        required=True)

    parser.add_argument("--out_header", action='store_true',
                        help="Set to include a header for output table.",
                        required = False, default = False)

    args = parser.parse_args()

    seqs = read_fasta(args.input)


    if args.out_header:
        print("\t".join(['name', 'canonical_start_codon_missing', 'canonical_stop_codon_missing', 'start_position',
                         'leading_percent_truncated', 'premature_stop_position', 'expected_stop_position',
                         'trailing_percent_truncated']))       

    for seq_id, seq in seqs.items():

        seq = seq.replace('-', '')

        seq = seq.upper()

        (exp_start_codon_present, start_codon_position, leading_percent_truncated) = start_codon_present(seq)

        (final_stop_present, premature_stop_codon_position, expected_stop_pos, trailing_percent_truncated) = stop_codon_premature_present(seq, start_codon_position)

        if not exp_start_codon_present or not final_stop_present or not math.isnan(premature_stop_codon_position):

            print("\t".join([args.name, str(not exp_start_codon_present), str(not final_stop_present), str(start_codon_position),
                             str(leading_percent_truncated), str(premature_stop_codon_position), str(expected_stop_pos), str(trailing_percent_truncated)]))


if __name__ == '__main__':

    main()

