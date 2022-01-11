#!/usr/bin/python3

import argparse
import os
import sys
from functions.io_utils import read_fasta

def main():

    parser = argparse.ArgumentParser(

            description="Parse coordinates from header lines of "
                        "prodigal-formatted FASTA. Return bedfile, where "
                        "optionally the coordinates are trimmed by a certain "
                        "amount. Can adjust min size of trimmed genes to keep.",

    epilog='''Usage example:

    python prodigal_fasta_to_trimmed_bed.py --input FASTA --trim_size 150 --min_size 100 --output BED

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="FASTA", type=str,
                        help="Path to input FASTA (with prodigal-formatted header lines)",
                        required=True)

    parser.add_argument("-o", "--output", metavar="BED", type=str,
                        help="Path to output BED", required=True)

    parser.add_argument("-t", "--trim_size", metavar="INT", type=int,
                        help="Number of bases to trim from the sides of each gene.",
                        required=False, default = 0)

    parser.add_argument("-m", "--min_size", metavar="INT", type=int,
                        help="Minimum size of genes to retain for bedfile.",
                        required=False, default = 1)

    parser.add_argument("-r", "--remove_trailing_scaffold_field",
                        action='store_true',
                        help="Set to remove last field delimited by \"_\" in scaffold name",
                        required = False, default = False)

    args = parser.parse_args()

    seqs = read_fasta(args.input, cut_header = False)

    with open(args.output, "w") as output_bed:
        for seqid in seqs.keys():
            seqid_split = seqid.split()

            scaffold = seqid_split[0].replace(">", "")
            start = int(seqid_split[2]) - 1
            stop = int(seqid_split[4])
            seq_name = seqid_split[8]

            start += args.trim_size
            stop -= args.trim_size

            final_size = stop - start

            if final_size >= args.min_size:

                if args.remove_trailing_scaffold_field:
                    scaffold = "_".join(scaffold.split("_")[0:-1])

                print(scaffold + "\t" + str(start) + "\t" + str(stop) + "\t" + seq_name,
                      file = output_bed)

            else:
                print("Dropping this line due to small size:\n" +
                      seqid + "\n\n", file = sys.stderr)


if __name__ == '__main__':
    main()
