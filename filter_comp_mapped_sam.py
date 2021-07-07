#!/usr/bin/python3

import argparse
import os
import sys

def main():

    parser = argparse.ArgumentParser(

            description="Filter SAM file that represents competitive mapping "
                        "to a database. Requires a FASTA file of the "
                        "reference sequences you are actually interested in "
                        "(i.e., those that are not just there to be "
                        "competitve mapping targets). Will parse SAM and only "
                        "keep reads mapped to the target sequences. Will also "
                        "filter out any reads that mapped to the target "
                        "sequences non-uniquely or where at least one of the "
                        "forward/reverse reads mapped to a non-target. In "
                        "addition to the filtered SAM, a breakdown of how "
                        "many reads multi-mapped to the target sequences, how often the forward and reverse 

    epilog='''Usage example:

    python add_Ns_around_seqs.py --input FASTA --output OUT_FASTA

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="FASTA", type=str,
                        help="Path to input FASTA", required=True)

    parser.add_argument("-o", "--output", metavar="FASTA", type=str,
                        help="Path to input FASTA", required=True)

    args = parser.parse_args()

    seqs = read_fasta(args.input, cut_header = False)

    seqs_padded = {}

    for seqid, sequence in seqs.items():

        seqs_padded[seqid] = 50 * "N" + sequence + 50 * "N"

    write_fasta(seqs_padded, args.output)


if __name__ == '__main__':

    main()
