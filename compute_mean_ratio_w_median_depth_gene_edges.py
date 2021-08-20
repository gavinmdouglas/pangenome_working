#!/usr/bin/python3

import argparse
import os
import sys
from statistics import median
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(

            description="This is a specialized script meant for getting a breakdown of how coverage drops off at the beginning and ending of scaffolds (which are individual genes in this case). "
                        "Specifically it parses per-site bedtools coverage files and outputs the average ratio of depth per site to the median per gene (based on all sites with non-zero depth only). "
                        "This is reported as the average across all genes in the first N and trailing N positions in each gene. "
                        "Importantly this script assumes that the reads were mapped to scaffolds where each scaffold was a different gene and so it makes sense that reads at the beginning and end are less commonly mapped. "
                        "It also assumes that the positions for each gene are in order (so that for a gene of 200 bp in length that the positions range from 1 to 200 in order in the input bedGraph. "
                        "Last any genes shorter than N * 2 bp or where < 80\\% of sites are covered by at least one read will be ignored. ",

    epilog='''Usage example:

    python compute_mean_ratio_w_median_depth_gene_edges.py --input per_site.bedGraph --output OUTPUT

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="per_site.bedGraph", type=str,
                        help="Path to input per-site bedGraph file created by bedtools coverage (with the -d option).",
                        required=True)

    parser.add_argument("-n", "--number_of_sites", metavar="INT", type=int,
                        help="Number of sites upsteam and downstream to consider for the output",
                        required=False, default = 250)

    parser.add_argument("-o", "--output", metavar="OUTPUT", type=str,
                        help="Path to output table.", required=True)

    parser.add_argument('--verbose', default=False, action='store_true',
                        help='If specified, print IDs of genes that are too short and will be ignored.')

    args = parser.parse_args()

    per_site_depth = defaultdict(list)

    with open(args.input, 'r') as BEDGRAPH:

        for line in BEDGRAPH:

            line = line.rstrip()

            line_split = line.split()

            gene_name = line_split[0]
            site_depth = int(line_split[4])

            # Since bp should be in order
            per_site_depth[gene_name].append(site_depth)


    summed_depth_ratio = [0] * args.number_of_sites * 2

    num_genes = 0

    for gene_name in per_site_depth.keys():

        if len(per_site_depth[gene_name]) < args.number_of_sites * 2:

            if args.verbose:
                print("Gene " + gene_name + " skipped because it is shorter than the (specified number of sites) * 2.",
                      file = sys.stderr)
            continue

        num_nonzero_pos = 0
        nonzero_pos_depth = []

        for pos in range(len(per_site_depth[gene_name])):
            if per_site_depth[gene_name][pos] > 0:
                num_nonzero_pos += 1
                nonzero_pos_depth.append(per_site_depth[gene_name][pos])

        if num_nonzero_pos / len(per_site_depth[gene_name]) < 0.8:
            continue

        num_genes += 1

        gene_median = median(nonzero_pos_depth)

        for i in range(args.number_of_sites):

            i_neg = (i + 1) * -1

            summed_depth_ratio[i] += per_site_depth[gene_name][i] / gene_median
            summed_depth_ratio[i_neg] += per_site_depth[gene_name][i_neg] / gene_median

    mean_depth_ratio = [x / num_genes for x in summed_depth_ratio]

    with open(args.output, 'w') as OUTPUT:

        print("gene_position\tmean_ratio", file = OUTPUT)

        for i in range(args.number_of_sites):
            i_1_based = str(i + 1)
            print(i_1_based + "\t" + str(mean_depth_ratio[i]), file = OUTPUT)

        for i in range(args.number_of_sites):
            i_neg = (i + 1) * -1
            print(str(i_neg) + "\t" + str(mean_depth_ratio[i_neg]), file = OUTPUT)


if __name__ == '__main__':
    main()
