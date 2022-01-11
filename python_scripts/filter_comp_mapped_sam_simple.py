#!/usr/bin/python3

import argparse

def main():

    parser = argparse.ArgumentParser(

            description="Filter SAM file that represents competitive mapping "
                        "to a database. Requires a FASTA file of the "
                        "reference sequences you are actually interested in "
                        "(i.e., those that are not just there to be "
                        "competitve mapping targets). Will parse SAM and only "
                        "keep reads mapped to the target sequences. Note that "
                        "this differs from the earlier version of the script "
                        "(non-simple) that did extensive purging of sequences "
                        "with any secondary sequences. This was will likely "
                        "keep more false alignments while keeping more reads "
                        "in general. All header lines are kept.",

    epilog='''Usage example:

    python filter_comp_mapped_sam_simple.py --input SAM --fasta FASTA --output SAM

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="SAM", type=str,
                        help="Path to input SAM file", required=True)

    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        help="Path to input FASTA file with target sequences "
                             "of interest (i.e., that contains header-lines "
                             "that correspond to the actual reference "
                             "sequences you are interested in getting "
                             "reliable read mappings to.", required=True)

    parser.add_argument("-o", "--output", metavar="SAM", type=str,
                        help="Path to output (filtered) SAM file",
                        required=True)

    args = parser.parse_args()

    # Parse FASTA and get target sequences.
    target_seqs = set()
    with open(args.fasta, "r") as fasta_in:
        for fasta_line in fasta_in:
            if fasta_line[0] == ">":
                target_seqs.add(fasta_line[1:].rstrip())
 
    sam_output = open(args.output, "w")

    with open(args.input, "r") as sam_read:
        for read_line in sam_read:
            read_line_split = read_line.split()
            
            # Print out headers.
            if read_line_split[0] in ["@HD", "@PG", "@SQ", "@RG", "@CO"]:
                print(read_line, end = "", file = sam_output)
                continue

            ref_hit = read_line_split[2]
            if ref_hit in target_seqs:
                print(read_line, end = "", file = sam_output)

    sam_output.close()


if __name__ == '__main__':
    main()
