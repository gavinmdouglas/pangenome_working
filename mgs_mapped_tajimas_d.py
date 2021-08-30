#!/usr/bin/python3

import argparse
import os
import sys
import pysam
import scipy.special
from statistics import mean, median
from math import ceil
from random import sample
from functions.pop_gen import (num_pairwise_diff_bases,
                               tajimas_d_and_wattersons_theta)

from functions.io_utils import read_vcf_variant_bases


def main():

    parser = argparse.ArgumentParser(

            description="Compute Tajima's D, Watterson's theta, and nucleotide diversity based on reads aligned to a reference genome.",

    epilog='''Usage example:

    python mgs_mapped_tajimas_d.py --bam BAM  --bed BED --vcf VCF --output OUTPUT_TABLE

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--bam", metavar="BAM", type=str, required=True,
                        help="Path to input BAM file.")

    parser.add_argument("--bed", metavar="BED", type=str, required=True,
                        help="BED containing coordinates of genes in reference "
                             "that was mapped to in the input BAM. Each line "
                             "should be a different gene and must include a "
                             "gene name.")

    parser.add_argument("--vcf", metavar="VCF", type=str, required=True,
                        help="Single-sample VCF containing the base call at "
                             "variant sites that should be used. This is "
                             "important because reads with bases at these "
                             "sites wont be counted at this position due to "
                             "the risk of them being sequencing errors (they "
                             "will just be ignored).")

    parser.add_argument("--rare_depth", metavar="INT", type=int,
                        required=False, default = None,
                        help="If specified, randomly subsample reads at each "
                             "site to specified depth before identifying "
                             "variants and differences.")

    parser.add_argument("-o", "--output", metavar="OUTPUT", type=str,
                        help="Path to output table.", required=True)

    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam, 'rb')

    outhandle = open(args.output, 'w')

    print("\t".join(["gene", "contig", "start_0_in", "stop_0_ex",
                     "mean_coverage", "median_coverage", "num_mapped_pos",
                     "breadth_coverage_nonzero", "breadth_coverage_atleast2",
                     "num_pairwise_diff", "num_comparisons",
                     "num_segregating_sites", "mean_polymorphic_coverage",
                     "theta_pi", "wattersons_theta", "tajimas_d"]),
                     file = outhandle)

    observed_variable_bases = read_vcf_variant_bases(in_vcf = args.vcf,
                                                     only_poly = True)

    if args.rare_depth:
        min_coverage = args.rare_depth
    else:
        min_coverage = 2

    with open(args.bed, 'r') as bedfile:

        for line in bedfile:
            
            line = line.rstrip()

            line_split = line.split()

            contig_name = line_split[0]
            gene_start_coor = int(line_split[1])
            gene_end_coor = int(line_split[2])
            gene_name = line_split[3]

            num_segregating_sites = 0
            num_comparisons = 0
            num_pairwise_diff = 0
            per_pos_coverage = list()
            poly_pos_coverage = list()
            num_nonzero_bases = 0
            num_atleastminread_bases = 0

            for pileupcolumn in bam.pileup(contig_name, gene_start_coor, gene_end_coor, stepper = 'nofilter', truncate = True):

                # Skip any positions not in the gene range.
                if pileupcolumn.pos < gene_start_coor or pileupcolumn.pos > gene_end_coor:
                    continue

                position_coor = contig_name + "|" + str(pileupcolumn.pos)

                bases_at_variable_site = []

                if position_coor in observed_variable_bases:
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            # query position is None if is_del or is_refskip is set.

                            query_base = pileupread.alignment.query_sequence[pileupread.query_position].upper()

                            # Only count mappings that match a base reported in the VCF.
                            if query_base in observed_variable_bases[position_coor]:
                                bases_at_variable_site.append(query_base)

                    pos_coverage = len(bases_at_variable_site)

                    if args.rare_depth and pos_coverage >= min_coverage:
                        bases_at_variable_site = sample(bases_at_variable_site,
                                                        min_coverage)
                        pos_coverage = min_coverage

                else:

                    if args.rare_depth:
                        pos_coverage = min_coverage
                    else:
                        pos_coverage = len(pileupcolumn.pileups)

                if pos_coverage == 0:
                    if position_coor in observed_variable_bases:
                        sys.exit("Error - site " + position_coor +
                                 " has 0 coverage but was expected to be a variant based on VCF.")
                else:
                    num_nonzero_bases += 1

                    if pos_coverage >= min_coverage:
                        num_atleastminread_bases += 1
                        num_comparisons += scipy.special.comb(pos_coverage, 2)

                        if bases_at_variable_site and len(set(bases_at_variable_site)) > 1:
                            pos_pairwise_diff = num_pairwise_diff_bases(bases_at_variable_site)
                            num_pairwise_diff += pos_pairwise_diff

                            if pos_pairwise_diff == 0:
                                sys.exit("Error - site " + position_coor +
                                         " has no pairwise differences but was expected to be a variant based on VCF.")

                            num_segregating_sites += 1
                            poly_pos_coverage.append(pos_coverage)

                per_pos_coverage.append(pos_coverage)

            if len(per_pos_coverage) > 0:
                mean_coverage = mean(per_pos_coverage)
                median_coverage = median(per_pos_coverage)
                nonzero_breadth = (num_nonzero_bases / (gene_end_coor - gene_start_coor)) * 100
                min_cov_breadth = (num_atleastminread_bases / (gene_end_coor - gene_start_coor)) * 100
            else:
                mean_coverage = 0
                median_coverage = 0
                nonzero_breadth = 0
                min_cov_breadth = 0

            if len(poly_pos_coverage) > 0:
                rounded_mean_poly_coverage = ceil(mean(poly_pos_coverage))
                nulc_div = (num_pairwise_diff / num_comparisons) * num_atleastminread_bases

                if rounded_mean_poly_coverage > 3:
                    (tajimas_d, wattersons_theta) = tajimas_d_and_wattersons_theta(theta_pi=nulc_div,
                                                                                   num_seg_sites=num_segregating_sites,
                                                                                   n=rounded_mean_poly_coverage)
                else:
                    tajimas_d = "NA"
                    wattersons_theta = num_segregating_sites / sum(1 / x for x in range(1, rounded_mean_poly_coverage))
            else:
                rounded_mean_poly_coverage = "NA"
                nulc_div = "NA"
                tajimas_d = "NA"
                wattersons_theta = "NA"

            outline = [gene_name, contig_name, gene_start_coor, gene_end_coor,
                       mean_coverage, median_coverage, num_nonzero_bases,
                       nonzero_breadth, min_cov_breadth, num_pairwise_diff,
                       num_comparisons, num_segregating_sites,
                       rounded_mean_poly_coverage, nulc_div, wattersons_theta,
                       tajimas_d]

            outline = [str(x) for x in outline]

            print("\t".join(outline), file = outhandle)

    outhandle.close()

if __name__ == '__main__':
    main()
