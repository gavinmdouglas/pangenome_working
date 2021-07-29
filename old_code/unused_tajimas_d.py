### Working commands not used because I realized just using the mpileup information would be easier than actually looking at individual reads ###

# for i in range(1):

#         mapped_reads = []

#         # Clean up and get positional and sequence info for all reads that mapped to a gene.
#         for read in sam.fetch(contig_name, gene_start_coor, gene_end_coor):

#             read_info = str(read).split()
#             map_start = int(read_info[3])
#             read_seq = read_info[9]
#             map_end = map_start + len(read_seq) - 1

#             mapped_read = aligned_info(contig_name, map_start, map_end, read_seq)

#             mapped_reads.append(trim_aligned_seq_to_element(mapped_read, gene_start_coor, gene_end_coor))

#         tajimas_d_info = read_based_tajimas_d_stats(mapped_reads)

#         print(tajimas_d_info)


# class aligned_info:
#     def __init__(self, contig_name, map_start, map_stop, seq):
#         self.contig_name = contig_name
#         self.map_start = map_start
#         self.map_stop = map_stop
#         self.seq = seq


# def trim_aligned_seq_to_element(mapped_read_info, element_start, element_stop):
#     '''Keep only subsequence of read/sequence that aligns to element. Assumes
#     that the positional coordinates refer to the same contig/chromosome.'''

#     if mapped_read_info.map_start < element_start:
#         leading_nonintersect_length = element_start - mapped_read_info.map_start
#         mapped_read_info.map_start = mapped_read_info.map_start + leading_nonintersect_length
#         mapped_read_info.seq = mapped_read_info.seq[leading_nonintersect_length:]

#     if mapped_read_info.map_stop > element_stop:
#         trailing_nonintersect_length = element_stop - mapped_read_info.map_stop
#         mapped_read_info.seq = mapped_read_info.seq[0:trailing_nonintersect_length]
#         mapped_read_info.map_stop = element_stop

#     return(mapped_read_info)


# def get_pairwise_seq_intersection(align_seq_a, align_seq_b):
#     '''Given two aligned_info objects, figure out if they intersect at all and
#     if so return a list of the portions of each sequence that do intersect.'''

#     # Quick check to see if they aligned to the same contig.
#     if align_seq_a.contig_name != align_seq_b.contig_name:
#         return(None)

#     if align_seq_a.map_start >= align_seq_b.map_start and align_seq_a.map_start <= align_seq_b.map_stop:
#         leading_diff = align_seq_a.map_start - align_seq_b.map_start
        
#         if align_seq_a.map_stop == align_seq_b.map_stop:
#             return([align_seq_a.seq,
#                     align_seq_b.seq[leading_diff:]])

#         elif align_seq_a.map_stop > align_seq_b.map_stop:
#             trailing_diff = align_seq_b.map_stop - align_seq_a.map_stop

#             return([align_seq_a.seq[:trailing_diff],
#                     align_seq_b.seq[leading_diff:]])

#         elif align_seq_a.map_stop < align_seq_b.map_stop:
#             trailing_diff = align_seq_a.map_stop - align_seq_b.map_stop

#             return([align_seq_a.seq,
#                     align_seq_b.seq[leading_diff:trailing_diff]])


#     elif align_seq_b.map_start >= align_seq_a.map_start and align_seq_b.map_start <= align_seq_a.map_stop:
#         leading_diff = align_seq_b.map_start - align_seq_a.map_start
        
#         if align_seq_a.map_stop == align_seq_b.map_stop:
#             return([align_seq_a.seq[leading_diff:],
#                     align_seq_b.seq])

#         elif align_seq_a.map_stop > align_seq_b.map_stop:
#             trailing_diff = align_seq_b.map_stop - align_seq_a.map_stop

#             return([align_seq_a.seq[leading_diff:trailing_diff],
#                     align_seq_b.seq])

#         elif align_seq_a.map_stop < align_seq_b.map_stop:
#             trailing_diff = align_seq_a.map_stop - align_seq_b.map_stop

#             return([align_seq_a.seq[leading_diff:],
#                     align_seq_b.seq[:trailing_diff]])
#     else:
#         return(None)


# def read_based_nucl_diversity(aligned_seqs):
#     '''Given a list of aligned_info objects, will determine what the total 
#     number of pairwise differences is. Will also return the the total number
#     of pairwise sites compared and the number of differences divided by this
#     quantity, which is the nucleotide diversity.'''

#     total_diff = 0
#     total_compare = 0

#     for i in range(len(aligned_seqs) - 1):
#         j = i + 1
#         while j <= len(aligned_seqs) - 1:

#             align_intersect = get_pairwise_seq_intersection(align_seq_a = aligned_seqs[i],
#                                                             align_seq_b = aligned_seqs[j])

#             if align_intersect is not None:

#                 if len(align_intersect[0]) != len(align_intersect[1]):
#                     sys.exit("Error - intersecting subsequences of aligned reads are of different sizes!")

#                 total_diff += sum(c1 != c2 for c1, c2 in zip(align_intersect[0],
#                                                              align_intersect[1]))

#                 total_compare += len(align_intersect[0])

#             j += 1

#     nucl_div = total_diff / total_compare

#     return((total_diff, total_compare, nucl_div))
