#!/usr/bin/python3

import argparse
import os
import sys

def main():

    parser = argparse.ArgumentParser(

            description="Convert inStrain profile SNV table to single-sample VCF",

    epilog='''Usage example:

    python instrain_profile_to_vcf.py --input PROFILE_SNVs --output VCF --default 5

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="PROFILE_SNVs", type=str,
                        help="Path to SNV table output by inStrain profile. command",
                        required=True)

    parser.add_argument("-m", "--min_depth", metavar="DEPTH", type=int,
                        help="Min depth to include an alternative allele in the VCF.",
                        required=False, default=1)

    parser.add_argument("-o", "--output", metavar="VCF", type=str,
                        help="Path to output VCF.", required=True)

    args = parser.parse_args()

    sample_name = os.path.splitext(os.path.basename(args.input))[0].replace("_SNVs", "")

    output_vcf = open(args.output, 'w')

    print("##fileformat=VCFv4.2", file = output_vcf)
    
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
          file = output_vcf)
    
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name,
          file = output_vcf)

    lc = 0
    with open(args.input, 'r') as input_SNVs:
        for line in input_SNVs:
            if lc == 0:
                lc += 1
                continue

            line_split = line.split()

            if int(line_split[2]) < args.min_depth or line_split[-1] == "AmbiguousReference":
                continue

            CHROM = line_split[0]
            POS = str(line_split[1] + 1)
            ID = CHROM + "|" + POS
            REF = line_split[4]

            bases = ['A', 'C', 'T', 'G']
            base_indices = [10, 11, 12, 13]
            alt_bases = []

            for i in range(4):
                base = bases[i]
                if base == REF:
                    ref_freq = int(line_split[base_indices[i]])
                    continue
                elif int(line_split[base_indices[i]]) >= args.min_depth:
                    alt_bases.append(base)

            ALT = ",".join(alt_bases)

            if len(alt_bases) == 0:
                continue

            GENOTYPE = []

            if ref_freq > args.min_depth:
                GENOTYPE.append("0")

            for j in range(1, len(alt_bases) + 1):                
                GENOTYPE.append(str(j))

            GENOTYPE = "/".join(GENOTYPE)

            print("\t".join([CHROM, POS, ID, REF, ALT, ".", ".", ".", "GT", GENOTYPE]),
                  file = output_vcf)

    output_vcf.close()

if __name__ == '__main__':
    main()
