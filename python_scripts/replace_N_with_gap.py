#!/usr/bin/python3

import argparse
import os
import sys
from functions.io_utils import read_fasta, write_fasta

def main():

    parser = argparse.ArgumentParser(

            description="Replace N characters with gaps in FASTA file.",

    epilog='''Usage example:

    python replace_N_with_gap.py -i IN.fasta -o OUT.fasta

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="FASTA", type=str,
                        help="Path to input FASTA.", required=True)

    parser.add_argument("--skip_if_no_Ns", action='store_true',
                        help="Set to only create output file if there is at least one N in original FASTA)",
                        required = False, default = False)

    parser.add_argument("-o", "--output", metavar="FASTA", type=str,
                        help="Path to output FASTA.", required=True)

    args = parser.parse_args()

    input_fasta = read_fasta(args.input)

    N_flag = False

    for seq_id in input_fasta.keys():

        if "N" in input_fasta[seq_id]:
            N_flag = True
            input_fasta[seq_id] = input_fasta[seq_id].replace("N", "-")

    if not args.skip_if_no_Ns or N_flag:
        write_fasta(input_fasta, args.output)


if __name__ == '__main__':
    main()
