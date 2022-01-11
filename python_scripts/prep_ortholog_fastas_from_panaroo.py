#!/usr/bin/python3

import argparse
import os
import sys
import textwrap
import gzip

def main():

    parser = argparse.ArgumentParser(

            description="Create FASTA files containing all genes classified by panaroo as part of each ortholog.",

    epilog='''Usage example:

    python prep_ortholog_fastas_from_panaroo.py -g gene_presence_absence.csv -d gene_data.csv -o OUTPUT_FOLDER

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-g", "--gene_presence_absence", metavar="TABLE", type=str,
                        help="Path to input gene presence absence table (gzipped) output by panaroo",
                        required=True)

    parser.add_argument("-d", "--gene_data", metavar="TABLE", type=str,
                        help="Path to input gene data table (gzipped) output by panaroo",
                        required=True)

    parser.add_argument("-p", "--prefix", metavar="STRING", type=str,
                        help="Prefix to add to output FASTA file names.",
                        required=True)

    parser.add_argument("-o", "--output", metavar="FOLDER", type=str,
                        help="Path to output folder, which will contain a FASTA for each ortholog listed in input gene presence absence table.",
                        required=True)

    args = parser.parse_args()

    ortholog_gene_sequences = dict()
    lc = 0
    with gzip.open(args.gene_data, 'rt') as in_gene_data:
        for line in in_gene_data:
            if lc == 0:
                lc += 1
                continue
            line = line.rstrip().split(',')
            ortholog_gene_sequences[line[3]] = line[5]

    os.mkdir(args.output)

    lc = 0
    with gzip.open(args.gene_presence_absence, 'rt') as in_presence_absence:
        for line in in_presence_absence:
            if lc == 0:
                lc += 1
                continue
            line = line.rstrip().split(',')

            with open(args.output + '/' + args.prefix + '_' + line[0] + '.fa', 'w') as ortholog_fasta:
                for gene in line[3:]:

                    if gene == '':
                        continue

                    for gene_split in gene.split(";"):

                        if '_stop' == gene_split[-5:]:
                            gene_split = gene_split.replace('_stop', '')

                        if '_len' == gene_split[-4:]:
                            gene_split = gene_split.replace('_len', '')

                        print(">" + gene_split, file = ortholog_fasta)
                        print(textwrap.fill(ortholog_gene_sequences[gene_split], width=70),
                              file = ortholog_fasta)

if __name__ == '__main__':
    main()
