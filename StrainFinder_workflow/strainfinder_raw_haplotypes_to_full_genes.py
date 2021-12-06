#!/usr/bin/python3

import argparse
import os
import sys
import gzip
import textwrap

def read_fasta(filename, cut_header=False):
    '''Read in FASTA file (gzipped or not) and return dictionary with each
    independent sequence id as a key and the corresponding sequence string as
    the value.'''

    # Intitialize empty dict.
    seq = {}

    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    # Read in FASTA line-by-line.
    if filename[-3:] == ".gz":
        fasta_in = gzip.open(filename, "rt")
    else:
        fasta_in = open(filename, "r")

    for line in fasta_in:

        line = line.rstrip()

        if len(line) == 0:
            continue

        # If header-line then split by whitespace, take the first element,
        # and define the sequence name as everything after the ">".
        if line[0] == ">":

            if cut_header:
                name = line.split()[0][1:]
            else:
                name = line[1:]

            name = name.rstrip("\r\n")

            # Make sure that sequence id is not already in dictionary.
            if name in seq:
                sys.stderr("Stopping due to duplicated id in file: " + name)

            # Intitialize empty sequence with this id.
            seq[name] = ""

        else:
            # Remove line terminator/newline characters.
            line = line.rstrip("\r\n")

            # Add sequence to dictionary.
            seq[name] += line

    fasta_in.close()

    return seq


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

    out_fasta = open(args.output, 'w')

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

                if h[variant_num] in struct_variants[var_pos].keys() and struct_variants[var_pos][h[variant_num]] != ref_allele:

                    obs_ref_allele = full_seq[site_i:site_i + len(ref_allele)]

                    if ref_allele != obs_ref_allele:
                        print("\n" + ref_allele + "\t" + obs_ref_allele, file = sys.stderr)
                        print("Structural ref alleles do not match at position " + str(var_pos), file = sys.stderr)
                        print("Current strategy is just to not add in this struc variant as it must be discordant with another variant.", file = sys.stderr)
                        #sys.exit("Structural ref alleles do not match at position " + str(var_pos))
                    else:
                        full_seq = full_seq[0:site_i] + struct_variants[var_pos][h[variant_num]] + full_seq[site_i + len(ref_allele):]

            else:

                full_seq = full_seq[0:site_i] + h[variant_num] + full_seq[site_i + 1:]

            variant_num -= 1

        out_fasta.write(">haplotype" + str(haplotype_num) + "_" + args.gene_name.replace(".fa", "") + "\n")
        out_fasta.write(textwrap.fill(full_seq, width=70) + "\n")

        haplotype_num += 1

    out_fasta.close()


if __name__ == '__main__':
    main()
