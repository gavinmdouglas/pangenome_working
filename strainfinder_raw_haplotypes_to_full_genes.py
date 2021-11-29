#!/usr/bin/python3

import argparse
import os
import sys
from functions.io_utils import read_fasta, write_fasta

def main():

    parser = argparse.ArgumentParser(

            description="Convert haplotypes in StrainFinder OTU tables to FASTA with invariant sites included and structural variants swapped in for placeholder nucleotides (if specified).",

    epilog='''Usage example:

    python strainfinder_raw_haplotypes_to_full_genes.py -g GENE_NAME --otutable STRAINFINDER_OTUTABLE --sites SITES_FILE --struct_map STRUCT_MAP -r REF_FASTA -o OUT_FASTA

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-g", "--gene_name", metavar="GENE_NAME", type=str,
                        help="Gene name in FASTA to use as reference",
                        required=True)

    parser.add_argument("-r", "--ref_fasta", metavar="REF_FASTA", type=str,
                        help="Path to FASTA that contains representative genes.",
                        required=True)

    parser.add_argument("--otutable", metavar="TABLE", type=str,
                        help="Path to StrainFinder OTU table that includes the haplotypes as the header",
                        required=True)

    parser.add_argument("--sites", metavar="SITES_FILE", type=str,
                        help="Path to sites file that represents all variant sites in StrainFinder output.",
                        required=True)

    parser.add_argument("--struct_map", metavar="STRUCT_MAP", type=str,
                        help="Path to file that represents the structural variants and mappings of nucleotides to structural alleles that should be swapped in.",
                        required=False, default=None)

    parser.add_argument("-o", "--output", metavar="OUT_FASTA", type=str,
                        help="Path to output FASTA.",
                        required=True)

    args = parser.parse_args()

    rep_genes = read_fasta(args.ref_fasta)
    original_gene = rep_genes[args.gene_name].upper()
    del rep_genes

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

    if args.struct_map:
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


    haplotype_seqs = dict()
    haplotype_num = 1

    for h in haplotypes:
        
        full_seq = original_gene

        # Note that as the structural variants can change the overall sequence
        # size, it's important to move from the end of the sequence to the
        # beginning so the indices don't get screwed up.
        variant_num = len(haplotypes[0]) - 1

        for var_pos in sorted(variant_positions, reverse = True):
            
            site_i = var_pos - 1

            if var_pos in struct_variants:

                ref_allele = struct_variant_refs[var_pos]

                obs_ref_allele = full_seq[site_i:site_i + len(ref_allele)]

                if ref_allele != obs_ref_allele:
                    print(ref_allele + "\t" + obs_ref_allele)
                    sys.exit("Structural ref alleles do not match at position " + str(var_pos))

                full_seq = full_seq[0:site_i] + struct_variants[var_pos][h[variant_num]] + full_seq[site_i + len(ref_allele):]

            else:

                full_seq = full_seq[0:site_i] + h[variant_num] + full_seq[site_i + 1:]

            variant_num -= 1

        haplotype_seqs["haplotype" + str(haplotype_num) + "_" + args.gene_name.replace(".fa", "")] = full_seq

        haplotype_num += 1

    write_fasta(haplotype_seqs, args.output)


if __name__ == '__main__':
    main()
