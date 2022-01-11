#!/usr/bin/python3

import argparse
import re
import os
import sys
import textwrap

from collections import defaultdict

from functions.io_utils import read_fasta

def main():

    parser = argparse.ArgumentParser(

            description="Format gene FASTA like prodigal FASTA (assumes that each gene is on a different scaffold, which it will use as the name)",

    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        help="Path to gene FASTA", required=True)

    parser.add_argument("-o", "--output", metavar="OUTPUT", type=str,
                        help="Path to output FASTA file, with headers formatted as in Prodigal output (except for the ID field, which is just the gene name)",
                        required=True)

    args = parser.parse_args()

    seqs = read_fasta(args.fasta, cut_header = True)

    scaffold_count = defaultdict(int)

    with open(args.output, "w") as OUTPUT:

        for seq_id, sequence in seqs.items():

            print(">" + seq_id + " # 1 # " + str(len(sequence)) + " # 1 # " + seq_id, file = OUTPUT)

            print(sequence, file = OUTPUT)


if __name__ == '__main__':

    main()
