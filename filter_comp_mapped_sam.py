#!/usr/bin/python3

import argparse
import os
import sys
from collections import defaultdict


# Taken from https://www.biostars.org/p/13051/
def asbin(n):
    '''converted a number to its binary representative (padded with 0s)'''
    return str(bin(n))[2:].zfill(17)


# class human_readable_sam_flags:
#     def __init__(self, flag):
#         self.read_paired = bool(int(asbin(flag)[-1]))
#         self.read_mapped_in_proper_pair = bool(int(asbin(flag)[-2]))
#         self.read_unmapped = bool(int(asbin(flag)[-3]))
#         self.mate_unmapped = bool(int(asbin(flag)[-4]))
#         self.read_reverse_strand = bool(int(asbin(flag)[-5]))
#         self.mate_reverse_strand = bool(int(asbin(flag)[-6]))
#         self.first_in_pair = bool(int(asbin(flag)[-7]))
#         self.second_in_pair = bool(int(asbin(flag)[-8]))
#         self.not_primary_alignment = bool(int(asbin(flag)[-9]))
#         self.read_fails_platform_vendor_quality_checks = bool(int(asbin(flag)[-10]))
#         self.read_is_PCR_or_optical_duplicate = bool(int(asbin(flag)[-11]))
#         self.supplementary_alignment = bool(int(asbin(flag)[-12]))


class pertinent_human_readable_sam_flags:
    def __init__(self, flag):
        self.read_paired = bool(int(asbin(flag)[-1]))
        self.read_unmapped = bool(int(asbin(flag)[-3]))
        self.mate_unmapped = bool(int(asbin(flag)[-4]))
        self.first_in_pair = bool(int(asbin(flag)[-7]))
        self.second_in_pair = bool(int(asbin(flag)[-8]))
        self.not_primary_alignment = bool(int(asbin(flag)[-9]))


def determine_secondary_match(info, primary_hit, focal_seqs, secondary_matches):
    '''Determine whether a read has a secondary match and if so return where
    the second match was (e.g., same target sequence, different subcluster,
    etc.).'''

    if not info in secondary_matches:
        return("NA")
    
    secondary_hit = secondary_matches[info]

    if not secondary_matches[info] in focal_seqs:
        return("nontarget")

    elif not primary_hit in focal_seqs and secondary_matches[info] in focal_seqs:
        return("target")

    elif primary_hit in focal_seqs and secondary_matches[info] in focal_seqs:

        if primary_hit == secondary_matches[info]:
            return("same.target")
        elif check_cluster_id(primary_hit, secondary_matches[info]):
            return("same.subcluster")
        else:
            return("diff.target")


def determine_paired_info(read_name, first_in_pair_flag, second_in_pair_flag,
                          primary_hit, focal_seqs, all_primary_hits,
                          secondary_matches):
    '''For a paired end read determine where the primary alignment of the other
    read occurred (if at all).'''

    read_pair_info = ",".join([read_name,
                               "True",
                               str(not first_in_pair_flag),
                               str(not second_in_pair_flag)])

    if not read_pair_info in all_primary_hits:
        return("unmapped")

    elif not all_primary_hits[read_pair_info] in focal_seqs:
        return("nontarget")

    elif all_primary_hits[read_pair_info] in focal_seqs and read_pair_info in secondary_matches:
        return("target.with.secondary")

    elif all_primary_hits[read_pair_info] in focal_seqs and not read_pair_info in secondary_matches:

        if not primary_hit in focal_seqs:
            return("target")
        elif primary_hit == all_primary_hits[read_pair_info]:
            return("same.target")
        elif check_cluster_id(primary_hit, all_primary_hits[read_pair_info]):
            return("same.subcluster")
        else:
            return("diff.target")


def check_cluster_id(match1, match2):
    '''Check whether two strings are identical after removing suffix that
    starts with "_cN", where N is an integer.'''
                                
    if "_c" in match1 and "_c" in match2:
        match1_trimmed = "_".join(match1.split("_")[0:-1])
        match2_trimmed = "_".join(match2.split("_")[0:-1])

        if match1_trimmed == match2_trimmed:
            return(True)
    return(False)
                         

def main():

    parser = argparse.ArgumentParser(

            description="Filter SAM file that represents competitive mapping "
                        "to a database. Requires a FASTA file of the "
                        "reference sequences you are actually interested in "
                        "(i.e., those that are not just there to be "
                        "competitve mapping targets). Will parse SAM and only "
                        "keep reads mapped to the target sequences. Will also "
                        "filter out any reads that mapped to the target "
                        "sequences non-uniquely or where at least one of the "
                        "forward/reverse reads mapped to a non-target. In "
                        "addition to the filtered SAM, a breakdown of how "
                        "many reads multi-mapped to the target sequences. "
                        "Will also provide a breakdown on how often the "
                        "forward and reverse reads map to the same target "
                        "sequence, another target sequence, or non-target "
                        "sequences. Note that this script assumes that each "
                        "read can be aligned a maximum of two time (e.g., "
                        "with the bowtie2 -k 2 option). Also, if target "
                        "sequences (i.e., focal reference sequences as "
                        "identified in the FASTA) are identical except for "
                        "different suffixes in of format _cN (where N is an "
                        "integer), then they are assumed to be different "
                        "sub-clusters of the same larger clusters.",

    epilog='''Usage example:

    python filter_comp_mapped_sam.py --input FASTA --output OUT_FASTA

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

    parser.add_argument("-r", "--report", metavar="OUTFILE", type=str,
                        help="Path to output report file.",
                        required=True)

    args = parser.parse_args()

    # Parse FASTA and get target sequences.
    target_seqs = set()
    with open(args.fasta, "r") as fasta_in:
        for fasta_line in fasta_in:
            if fasta_line[0] == ">":
                target_seqs.add(fasta_line[1:].rstrip())
 

    # Parse SAM file and keep track of for all alignments called as secondary
    # alignments which reference sequence they hit. Keep this information
    # separate for single/paired-end reads and also for the first and second
    # mates of paired-end reads.
    secondary_hits = {}

    # Also for primary paired-end alignments, also keep track of what the hit
    # location was.
    primary_hits = {}

    with open(args.input, "r") as sam_initial:
        for initial_line in sam_initial:
            initial_line_split = initial_line.split()
            
            # Skip header-lines for initial read-through.
            if initial_line_split[0] in ["@HD", "@PG", "@SQ"]:
                continue

            read_id = initial_line_split[0]
            map_info = pertinent_human_readable_sam_flags(int(initial_line_split[1]))
            ref_hit = initial_line_split[2]

            read_info = ",".join([read_id,
                                  str(map_info.read_paired),
                                  str(map_info.first_in_pair),
                                  str(map_info.second_in_pair)])

            if map_info.not_primary_alignment:
                secondary_hits[read_info] = ref_hit

            elif map_info.read_paired:
                primary_hits[read_info] = ref_hit


    # Read through SAM file again, but this time write out lines that should
    # be retained. Also, keep track of read counts of different categories
    # to write on to report table.
    sam_output = open(args.output, "w")

    category_counts = defaultdict(int)

    with open(args.input, "r") as sam_reread:
        for reread_line in sam_reread:
            reread_line_split = reread_line.split()
            
            # Print out headers except for sequence headerlines that do not
            # match target sequences.
            if reread_line_split[0] in ["@HD", "@PG"]:
                print(reread_line, end = "", file = sam_output)
                continue

            elif reread_line_split[0] == "@SQ":
                sequence_id = reread_line_split[1].replace("SN:", "")
    
                if sequence_id in target_seqs:
                    print(reread_line, end = "", file = sam_output)
                continue

            read_id = reread_line_split[0]
            map_info = pertinent_human_readable_sam_flags(int(reread_line_split[1]))
            ref_hit = reread_line_split[2]

            read_info = ",".join([read_id,
                                  str(map_info.read_paired),
                                  str(map_info.first_in_pair),
                                  str(map_info.second_in_pair)])

            # Only consider primary alignments on this read-through.
            if not map_info.not_primary_alignment:

                # Confirm that read was mapped (or skip it).
                if not map_info.read_unmapped:
                    category_counts["either|any|any"] += 1
                else:
                    continue

                if ref_hit in target_seqs:
                    on_target_flag = "target"
                else:
                    on_target_flag = "nontarget"

                secondary_map_match = determine_secondary_match(info = read_info,
                                                                primary_hit = ref_hit,
                                                                focal_seqs = target_seqs,
                                                                secondary_matches = secondary_hits)

                if not map_info.read_paired:
                    paired_info = "NA"

                else:
                    paired_info = determine_paired_info(read_name = read_id,
                                                        first_in_pair_flag = map_info.first_in_pair,
                                                        second_in_pair_flag = map_info.second_in_pair,
                                                        primary_hit = ref_hit,
                                                        focal_seqs = target_seqs,
                                                        all_primary_hits = primary_hits,
                                                        secondary_matches = secondary_hits)


                if on_target_flag == "target" and \
                   secondary_map_match == "NA" and \
                   paired_info in ["NA", "unmapped", "same.target",
                                   "same.subcluster", "diff.target"]:

                    print(reread_line, end = "", file = sam_output)

                count_category = "|".join([on_target_flag, secondary_map_match,
                                           paired_info])

                category_counts[count_category] += 1

    sam_output.close()


    # Write out breakdown of read mappings to report file.
    with open(args.report, "w") as report_out:

        print("primary_hit|secondary_hit|paired_hit" + "\t" + os.path.basename(args.input),
              file = report_out)

        print("either|any|any" + "\t" + str(category_counts["either|any|any"]),
              file = report_out)

        for target_setting in ['target', 'nontarget']:
            for secondary_map_setting in ['NA', 'nontarget', 'target',
                                          'same.target', 'same.subcluster',
                                          'diff.target']:
                for paired_map_setting in ['NA', 'unmapped', 'nontarget',
                                           'target.with.secondary', 'target',
                                           'same.target', 'same.subcluster',
                                           'diff.target']:

                    out_category = "|".join([target_setting, secondary_map_setting,
                                               paired_map_setting])

                    print(out_category + "\t" + str(category_counts[out_category]),
                          file = report_out)


if __name__ == '__main__':
    main()
