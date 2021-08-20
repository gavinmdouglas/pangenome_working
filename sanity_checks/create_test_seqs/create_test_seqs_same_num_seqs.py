#!/usr/bin/python3

from pprint import pprint
from random import sample, randint

bases = ['A', 'C', 'T', 'G']

starting_seqs = []

num_reps = 1000

sample_dists = [[1, 2, 3, 4, 5, 6, 7, 8, 9],
                [3, 4, 5, 6, 7],
                [1, 2, 8, 9]]

for rep in range(num_reps):

    starting_seq = ""

    for i in range(1000):
        rand_base = sample(bases, 1)[0]
        starting_seq += rand_base

    starting_seqs.append(starting_seq)

    final_seqs = [starting_seq] * 10

    num_subs = sample(list(range(1, 100)), 1)[0]

    sub_sites = sample(list(range(1000)), num_subs)

    sample_dist = sample(sample_dists, 1)[0]

    for sub_site in sub_sites:

        possible_subs = [x for x in bases if x != starting_seq[sub_site]]

        sub_base = sample(possible_subs, 1)[0]

        num_seqs_w_sub = sample(sample_dist, 1)[0]
        seqs_w_sub = sample(list(range(10)), num_seqs_w_sub)

        for seq_i in seqs_w_sub:
            final_seqs[seq_i] = final_seqs[seq_i][:sub_site] + sub_base + final_seqs[seq_i][sub_site + 1:]

    final_file = 'test_seqs/final_seqs/final_' + str(rep + 1) + '.fasta'
    with open(final_file, 'w') as final_out:

        for seq_i in range(10):
            print(">seq" + str(seq_i + 1), file = final_out)
            print(final_seqs[seq_i], file = final_out)

with open('test_seqs/starting_seqs.fasta', 'w') as starting_out:
    for rep in range(num_reps):
        print(">starting_seq_" + str(rep + 1), file = starting_out)
        print(starting_seqs[rep], file = starting_out)
