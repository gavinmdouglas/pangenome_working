#!/usr/bin/python3

from random import sample, randint
from math import floor
from numpy.random import choice

bases = ['A', 'C', 'T', 'G']

starting_seqs = []

num_reps = 1000

possible_gc_content = [0.25, 0.50, 0.75]
possible_num_seqs = [5, 10, 15, 20]
possible_num_subs = [5, 25, 50, 75, 100]

rep_metadata = open('test_seqs_for_tajimas_d_varied/rep_metadata.tsv', 'w')
print("rep\tnum_seqs\tgc_content\tmax_percent_seqs_w_sub\tS", file = rep_metadata)

for rep in range(num_reps):

    starting_seq = ""

    # Base proportions, determined by randomly sampled gc content.
    sampled_gc_content = sample(possible_gc_content, 1)[0]

    base_prop = [(1 - sampled_gc_content) / 2,
                 sampled_gc_content / 2,
                 (1 - sampled_gc_content) / 2,
                 sampled_gc_content / 2]

    # Starting sequences are 1000 bp in all cases.
    for i in range(1000):
        rand_base = choice(a = bases, p = base_prop)[0]
        starting_seq += rand_base

    starting_seqs.append(starting_seq)

    num_seqs = sample(possible_num_seqs, 1)[0]

    final_seqs = [starting_seq] * num_seqs

    num_subs = sample(possible_num_subs, 1)[0]

    sub_sites = sample(list(range(1000)), num_subs)

    abs_max_num_seqs = floor(num_seqs / 2)

    if abs_max_num_seqs > 1:
        max_num_seqs_w_sub = sample(list(range(1, floor(num_seqs / 2))), 1)[0]
    else:
        max_num_seqs_w_sub = 1

    max_percent_seqs_w_sub = (max_num_seqs_w_sub / num_seqs) * 100

    print("\t".join([str(rep + 1), str(num_seqs), str(sampled_gc_content),
                    "{:.2f}".format(max_percent_seqs_w_sub), str(num_subs)]),
          file = rep_metadata)

    for sub_site in sub_sites:

        possible_subs = []
        possible_sub_prob = []
        possible_sub_prob_sum = 0
        for i in range(4):
            if bases[i] != starting_seq[sub_site]:
                possible_subs.append(bases[i])
                possible_sub_prob.append(base_prop[i])
                possible_sub_prob_sum += base_prop[i]

        possible_sub_prob = [x / possible_sub_prob_sum for x in possible_sub_prob]

        sub_base = choice(a = possible_subs, p = possible_sub_prob)[0]

        num_seqs_w_sub = sample(list(range(1, max_num_seqs_w_sub + 1)), 1)[0]
        seqs_w_sub = sample(list(range(num_seqs)), num_seqs_w_sub)

        for seq_i in seqs_w_sub:
            final_seqs[seq_i] = final_seqs[seq_i][:sub_site] + sub_base + final_seqs[seq_i][sub_site + 1:]

    final_file = 'test_seqs_for_tajimas_d_varied/final_seqs/final_' + str(rep + 1) + '.fasta'
    with open(final_file, 'w') as final_out:

        for seq_i in range(num_seqs):
            print(">rep" + str(rep + 1) + "_seq" + str(seq_i + 1), file = final_out)
            print(final_seqs[seq_i], file = final_out)

with open('test_seqs_for_tajimas_d_varied/starting_seqs.fasta', 'w') as starting_out:
    for rep in range(num_reps):
        print(">starting_seq" + str(rep + 1), file = starting_out)
        print(starting_seqs[rep], file = starting_out)
