#!/usr/bin/python3

import argparse
import re
import os
import sys
import textwrap

from functions.io_utils import read_fasta, write_fasta

def main():

    parser = argparse.ArgumentParser(

            description="Add 50 Ns around all sequences in specified FASTA file.",

    epilog='''Usage example:

    python add_Ns_around_seqs.py --input FASTA --output OUT_FASTA

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="FASTA", type=str,
                        help="Path to input FASTA", required=True)

    parser.add_argument("-o", "--output", metavar="FASTA", type=str,
                        help="Path to output FASTA", required=True)

    args = parser.parse_args()

    seqs = read_fasta(args.input, cut_header = False)

    seqs_padded = {}

    for seqid, sequence in seqs.items():

        seqs_padded[seqid] = 50 * "N" + sequence + 50 * "N"

    write_fasta(seqs_padded, args.output)


if __name__ == '__main__':

    main()
