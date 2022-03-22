#!/usr/bin/python3

import argparse
import os
import sys
from functions.io_utils import read_fasta
from functions.codon import start_codon_present, stop_codon_premature_present
import pandas as pd
import numpy as np


def main():

    parser = argparse.ArgumentParser(

    description='Call canonical start / stop codon knockout mutations and premature stop '
                'codons in gene FASTAs. Assumes that each input sequence '
                '*should* start with start codon and end in stop codon.',

    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input', metavar='FASTA', type=str,
                        help='DNA sequence FASTA of input genes.',
                        required=True)

    parser.add_argument('-o', '--output', metavar='TABLE', type=str,
                        help='Output table to create.', required=True)

    args = parser.parse_args()

    seqs = read_fasta(args.input)

    seq_info = pd.DataFrame(columns = ['canonical_start_codon_missing', 'canonical_stop_codon_missing', 'start_position',
                                       'leading_percent_truncated', 'premature_stop_position', 'expected_stop_position', 'trailing_percent_truncated'],
                            index = seqs.keys())
       

    for seq_id, seq in seqs.items():

        seq = seq.replace('-', '')

        seq = seq.upper()

        (exp_start_codon_present, start_codon_position, leading_percent_truncated) = start_codon_present(seq)

        (final_stop_present, premature_stop_codon_position, expected_stop_pos, trailing_percent_truncated) = stop_codon_premature_present(seq, start_codon_position)

        seq_info.loc[seq_id, 'canonical_start_codon_missing'] = not exp_start_codon_present
        seq_info.loc[seq_id, 'canonical_stop_codon_missing'] = not final_stop_present
        seq_info.loc[seq_id, 'start_position'] = start_codon_position
        seq_info.loc[seq_id, 'leading_percent_truncated'] = leading_percent_truncated
        seq_info.loc[seq_id, 'premature_stop_position'] = premature_stop_codon_position
        seq_info.loc[seq_id, 'expected_stop_position'] = expected_stop_pos
        seq_info.loc[seq_id, 'trailing_percent_truncated'] = trailing_percent_truncated

    seq_info.to_csv(args.output, sep = "\t", header = False, na_rep='NA',
                    index_label = 'sequence')


if __name__ == '__main__':

    main()

