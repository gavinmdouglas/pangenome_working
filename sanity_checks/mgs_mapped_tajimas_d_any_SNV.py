#!/usr/bin/python3

import argparse
import os
import sys
import pysam
import scipy.special
from statistics import mean, median
from math import ceil
sys.path.insert(1, '/home/gdouglas/projects/honey_bee/scripts/pangenome_working/')

from functions.pop_gen import (num_pairwise_diff_bases,
                               tajimas_d_and_wattersons_theta)


def main():

    parser = argparse.ArgumentParser(

            description="Compute Tajima's D, Watterson's theta, and nucleotide diversity based on reads aligned to a reference genome. NOTE THAT THIS VERSION OF THE SCRIPT WILL CONSIDER ANY NON-REFERENCE BASE IN READS AS A SNV, WHEREAS I THINK CALLING SNVS BASED ON AT LEAST ~5 READS PROBABLY MAKES MUCH MORE SENSE TO LIMIT THE EFFECT OF SEQUENCING ERRORS ON THESE METRICS.",

    epilog='''Usage example:

    python mgs_mapped_tajimas_d.py --bam BAM  --bed BED --output OUTPUT_TABLE

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--bam", metavar="BAM", type=str, required=True,
                        help="Path to input BAM file.")

    parser.add_argument("--bed", metavar="BED", type=str, required=True,
                        help="BED containing coordinates of genes in reference "
                             "that was mapped to in the input BAM. Each line "
                             "should be a different gene and must include a "
                             "gene name.")

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
            num_atleast2read_bases = 0

            for pileupcolumn in bam.pileup(contig_name, gene_start_coor, gene_end_coor, stepper = 'nofilter', truncate = True):

                # Skip any positions not in the gene range.
                if pileupcolumn.pos < gene_start_coor or pileupcolumn.pos > gene_end_coor:
                    continue

                unique_mapped_bases = []

                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        # query position is None if is_del or is_refskip is set.

                        query_base = pileupread.alignment.query_sequence[pileupread.query_position].upper()

                        # Only count mappings that are unambiguous.
                        if query_base in ['A', 'C', 'T', 'G']:
                            unique_mapped_bases.append(query_base)

                pos_coverage = len(unique_mapped_bases)

                if pos_coverage > 0:
                    num_nonzero_bases += 1

                if pos_coverage > 1:
                    num_atleast2read_bases += 1
                    num_comparisons += scipy.special.comb(pos_coverage, 2)
                    pos_pairwise_diff = num_pairwise_diff_bases(unique_mapped_bases)
                    num_pairwise_diff += pos_pairwise_diff

                    if pos_pairwise_diff > 0:
                        num_segregating_sites += 1
                        poly_pos_coverage.append(pos_coverage)

                per_pos_coverage.append(pos_coverage)

            if len(per_pos_coverage) > 0:
                mean_coverage = mean(per_pos_coverage)
                median_coverage = median(per_pos_coverage)
                nonzero_breadth = (num_nonzero_bases / (gene_end_coor - gene_start_coor)) * 100
                min2_breadth = (num_atleast2read_bases / (gene_end_coor - gene_start_coor)) * 100
            else:
                mean_coverage = 0
                median_coverage = 0
                nonzero_breadth = 0
                min2_breadth = 0

            if len(poly_pos_coverage) > 0:
                rounded_mean_poly_coverage = ceil(mean(poly_pos_coverage))
                nulc_div = (num_pairwise_diff / num_comparisons) * num_atleast2read_bases

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
                       nonzero_breadth, min2_breadth, num_pairwise_diff,
                       num_comparisons, num_segregating_sites,
                       rounded_mean_poly_coverage, nulc_div, wattersons_theta,
                       tajimas_d]

            outline = [str(x) for x in outline]

            print("\t".join(outline), file = outhandle)

    outhandle.close()

if __name__ == '__main__':
    main()
