#!/usr/bin/python3

import argparse
import os

from functions.io_utils import read_fasta

def main():

    parser = argparse.ArgumentParser(

            description="Format gene FASTA to be in gene format for kpileup as described as necessary for StrainFinder.",

    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        help="Path to gene FASTA", required=True)

    parser.add_argument("-o", "--output", metavar="OUTPUT", type=str,
                        help="Path to directory where a single output gene-formatted file will be created for separate sequence in the input FASTA.",
                        required=True)

    args = parser.parse_args()

    seqs = read_fasta(args.fasta, cut_header = True)

    os.mkdir(args.output)

    for seq_id in seqs.keys():

        outfile = args.output + "/" + seq_id + ".gene"

        sequence = seqs[seq_id]

        with open(outfile, "w") as OUTPUT:

            print("\t".join([seq_id, seq_id, "1", str(len(sequence)), "+",
                             sequence]),
                  file = OUTPUT)

if __name__ == '__main__':

    main()
