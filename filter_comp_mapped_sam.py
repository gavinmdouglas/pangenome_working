#!/usr/bin/python3

import argparse
import os
import sys


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

    # Also for primary paired-end alignments, also keep track of what the hit
    # location was.

    singleend_secondary_hits = {}
    pairedend_firstmate_secondary_hits = {}
    pairedend_secondmate_secondary_hits = {}

    pairedend_firstmate_primary_hits = {}
    pairedend_secondmate_primary_hits = {}

    with open(args.input, "r") as sam_initial:
        for initial_line in sam_initial:
            initial_line_split = initial_line.split()
            
            # Skip header-lines for initial read-through.
            if initial_line_split[0] in ["@HD", "@PG", "@SQ"]:
                continue

            read_id = initial_line_split[0]
            map_info = pertinent_human_readable_sam_flags(int(initial_line_split[1]))
            ref_hit = initial_line_split[2]

            if map_info.not_primary_alignment:
                if not map_info.read_paired:
                    singleend_secondary_hits[read_id] = ref_hit

                else:
                    if map_info.first_in_pair:
                        pairedend_firstmate_secondary_hits[read_id] = ref_hit
                    elif map_info.second_in_pair:
                        pairedend_secondmate_secondary_hits[read_id] = ref_hit
                    else:
                        sys.exit("Error, paired-end read not marked as first "
                                 "or second in pair: " + initial_line)

            elif map_info.read_paired:
                if map_info.first_in_pair:
                    pairedend_firstmate_primary_hits[read_id] = ref_hit
                elif map_info.second_in_pair:
                    pairedend_secondmate_primary_hits[read_id] = ref_hit
                else:
                    sys.exit("Error, paired-end read not marked as first "
                                 "or second in pair: " + initial_line)

    # Read through SAM file again, but this time write out lines that should
    # be retained. Also, keep track of read counts of different categories
    # to write on to report table.

    sam_output = open(args.output, "w")

    total_unique_reads = 0

    singleend_unique_target_map = 0
    singleend_unique_offtarget_map = 0

    singleend_multi_exact_target_map = 0
    singleend_multi_diff_subcluster_map = 0
    singleend_multi_diff_target_map = 0
    singleend_multi_one_offtarget_map = 0
    singleend_multi_both_offtarget_map = 0

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

            # Only consider primary alignments on this read-through.
            if not map_info.not_primary_alignment:

                # Confirm that read was mapped (or skip it).
                if not map_info.read_unmapped:
                    total_unique_reads += 1
                else:
                    continue

                target_db_hit = ref_hit in target_seqs

                # For single-end read, write out if there is no secondary
                # alignment and the read mapped to a target sequence.
                if not map_info.read_paired:
                    
                    if not read_id in singleend_secondary_hits:

                        if target_db_hit:
                            print(reread_line, end = "", file = sam_output)
                            singleend_unique_target_map += 1
                        else:
                            singleend_unique_offtarget_map += 1

                    # If there was a secondary hit then figure out if it
                    # corresponded to exactly the same target sequence (if the
                    # current read mapped to the target), the same subcluster
                    # (after removing the _cN suffix and given that the current
                    # read mapped to the target), or whether one or both of
                    # them mapped offtarget.
                    else:

                        secondary_db_hit = singleend_secondary_hits[read_id] in target_seqs 

                        if target_db_hit:
                            if not secondary_db_hit:
                                singleend_multi_one_offtarget_map += 1
                            elif ref_hit == singleend_secondary_hits[read_id]:
                                singleend_multi_exact_target_map += 1
                            elif check_cluster_id(ref_hit, singleend_secondary_hits[read_id]):
                                singleend_multi_diff_subcluster_map += 1
                            else:
                                singleend_multi_diff_target_map += 1


                        elif secondary_db_hit:
                            singleend_multi_one_offtarget_map += 1
                        else:
                            singleend_multi_both_offtarget_map += 1

                # For paired-end read, need to do similar work!
                else:
                    pass

    sam_output.close()

    # Write out breakdown of read mappings to report file.
    with open(args.report, "w") as report_out:

        print("#category " + os.path.basename(args.input), file = report_out)
        print("total_unique_reads " + str(total_unique_reads), file = report_out)
        print("singleend_unique_target_map " + str(singleend_unique_target_map), file = report_out)
        print("singleend_unique_offtarget_map " + str(singleend_unique_offtarget_map), file = report_out)
        print("singleend_multi_exact_target_map " + str(singleend_multi_exact_target_map), file = report_out)
        print("singleend_multi_diff_subcluster_map " + str(singleend_multi_diff_subcluster_map), file = report_out)
        print("singleend_multi_diff_target_map " + str(singleend_multi_diff_target_map), file = report_out)
        print("singleend_multi_one_offtarget_map " + str(singleend_multi_one_offtarget_map), file = report_out)
        print("singleend_multi_both_offtarget_map " + str(singleend_multi_both_offtarget_map), file = report_out)


if __name__ == '__main__':
    main()
