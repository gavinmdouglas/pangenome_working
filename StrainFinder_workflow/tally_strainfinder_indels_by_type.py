#!/usr/bin/python3

import argparse
import os
import sys

def main():

    parser = argparse.ArgumentParser(

            description="Count INDELs identified in StrainFinder haplotype compared with the reference gene, broken down by frameshift/non-frameshift and deletion/insertion. Again, these inferences of deletion/insertion are all relative to the reference gene.",

    epilog='''Usage example:

    python tally_strainfinder_indels_by_type.py -g GENE_NAME --otutable STRAINFINDER_OTUTABLE --sites SITES_FILE --struct_map STRUCT_MAP > OUTPUT

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-g", "--gene_name", metavar="GENE_NAME", type=str,
                        help="Gene name to include in output table.",
                        required=True)

    parser.add_argument("--otutable", metavar="TABLE", type=str,
                        help="Path to StrainFinder OTU table that includes the haplotypes as the header",
                        required=True)

    parser.add_argument("--sites", metavar="SITES_FILE", type=str,
                        help="Path to sites file that represents all variant sites in StrainFinder output.",
                        required=True)

    parser.add_argument("--struct_map", metavar="STRUCT_MAP", type=str,
                        help="Path to file that represents the structural variants and mappings of nucleotides to structural alleles that should be swapped in.",
                        required=True)

    parser.add_argument("--out_header", action='store_true',
                        help="Set to include a header for output table.",
                        required = False, default = False)

    args = parser.parse_args()

    with open(args.otutable, 'r') as in_otutable:
        for otu_line in in_otutable:
            haplotypes = otu_line.split()
            break

    variant_positions = []
    with open(args.sites, 'r') as in_sites:
        for site_line in in_sites:
            variant_positions.append(int(site_line.split()[1]))

    struct_variants = dict()
    struct_variant_refs = dict()

    bases = ['A', 'C', 'G', 'T']

    with open(args.struct_map, 'r') as in_struct_map:
        
        for struct_line in in_struct_map:
            struct_line_split = struct_line.split()
            pos = int(struct_line_split[1])
            ref_allele = struct_line_split[2]
            variant_alleles = struct_line_split[3].split(',')

            struct_variants[pos] = dict()

            for i in range(len(variant_alleles)):
                struct_variants[pos][bases[i]] = variant_alleles[i]

            struct_variant_refs[pos] = ref_allele

    if args.out_header:
        print("gene\thaplotype\tnonframe_insert\tnonframe_del\tframe_insert\tframe_del")

    haplotype_seqs = dict()
    haplotype_num = 1

    for h in haplotypes:

        nonframe_del = 0
        nonframe_insert = 0
        frame_del = 0
        frame_insert = 0

        variant_num = 0

        for var_pos in variant_positions:
            
            site_i = var_pos - 1

            if var_pos in struct_variants:

                ref_allele = struct_variant_refs[var_pos]

                if not h[variant_num] in struct_variants[var_pos]:
                    print("Error - base not found in structured mapping file - base " + h[variant_num] + " at position " + str(var_pos) + ". Assuming that this is the reference base.",
                          file = sys.stderr)
                    alt_allele = ref_allele

                else:
                    alt_allele = struct_variants[var_pos][h[variant_num]]

                # Only an indel if they are not the same length.
                if len(ref_allele) != len(alt_allele):

                    abs_diff = abs(len(alt_allele) - len(ref_allele))

                    if len(alt_allele) > len(ref_allele):

                        if abs_diff % 3 != 0:
                            frame_insert += 1
                        else:
                            nonframe_insert += 1

                    else:
                        if abs_diff % 3 != 0:
                            frame_del += 1
                        else:
                            nonframe_del += 1
                 
            variant_num += 1


        if nonframe_insert > 0 or nonframe_del > 0 or frame_insert > 0 or frame_del > 0:
            outline = [args.gene_name, str(haplotype_num), str(nonframe_insert),
                       str(nonframe_del), str(frame_insert), str(frame_del)]

            print("\t".join(outline))

        haplotype_num += 1

if __name__ == '__main__':
    main()
