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

    parser.add_argument("-g", "--gene_ids", metavar="GENE_IDS", type=str,
                        help="Path to file with gene ids that should be retained.",
                        required=False, default = None)

    parser.add_argument("-o", "--output", metavar="OUTPUT", type=str,
                        help="Path to directory where a single output gene-formatted file will be created for separate sequence in the input FASTA.",
                        required=True)

    args = parser.parse_args()

    seqs = read_fasta(args.fasta, cut_header = True)

    if args.gene_ids is not None:
        genes_to_keep = set()
        with open(args.gene_ids, 'r') as GENE_IDS:
            for line in GENE_IDS:
                genes_to_keep.add(line.rstrip())

    os.mkdir(args.output)

    for seq_id in seqs.keys():

        outfile = args.output + "/" + seq_id + ".gene"

        # Skip if not in set of genes to retain (if option specified).
        if args.gene_ids is not None and seq_id not in genes_to_keep:
            continue

        sequence = seqs[seq_id]

        with open(outfile, "w") as OUTPUT:

            print("\t".join([seq_id, seq_id, "1", str(len(sequence)), "+",
                             sequence]),
                  file = OUTPUT)

if __name__ == '__main__':

    main()
