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

            description="Read in PROKKA gene FASTA and GFF files and output FASTA with header info in same format as Prodigal gene FASTAs (except the ID field will differ from the Prodigal default)",

    epilog='''Usage example:

    python prodigal_fasta_get_prokka_ids.py --fasta FASTA --gff GFF --output OUT_FASTA

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        help="Path to PROKKA gene file", required=True)

    parser.add_argument("-g", "--gff", metavar="GFF", type=str,
                        help="Path to PROKKA GFF3 file", required=True)

    parser.add_argument("-o", "--output", metavar="OUTPUT", type=str,
                        help="Path to output FASTA file, with headers formatted as in Prodigal output (except for the ID field)",
                        required=True)

    args = parser.parse_args()

    seqs = read_fasta(args.fasta, cut_header = True)

    scaffold_count = defaultdict(int)

    with open(args.output, "w") as OUTPUT:

        with open(args.gff, "r") as GFF:

            for line in GFF:

                # Stop reading in once reached the FASTA section of GFF.
                if line[0:7] == "##FASTA":
                    break

                # Skip all lines starting with comment characters.
                elif line[0:2] == "##":
                    continue

                else:
                    line_split = line.split()

                    seq_info = line_split[8]

                    seq_id = seq_info.split(";")[0][3:]

                    if seq_id in seqs.keys():

                        scaffold = line_split[0]
                        scaffold_count[scaffold] += 1
                        
                        line_begin = ">" + scaffold + "_" + str(scaffold_count[scaffold])
                        start = line_split[3]
                        stop = line_split[4]

                        if line_split[6] == "+":
                            strand = "1"
                        elif line_split[6] == "-":
                            strand = "-1"

                        print(" # ".join([line_begin, start, stop, strand, seq_info]),
                              file = OUTPUT)

                        print (textwrap.fill(seqs[seq_id], width=70),
                               file = OUTPUT)

if __name__ == '__main__':

    main()
