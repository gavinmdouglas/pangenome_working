#!/usr/bin/python3

from math import sqrt
import numpy as np
import itertools
import sys

def tajimas_d_and_diversity_metrics(in_seqs):
    '''Return n, S, Watterson's theta, nucleotide diversity, and Tajima's D as a tuple
    for a list of input sequences.'''

    (theta_pi, S) = calc_theta_pi_and_num_segregating(in_seqs)

    N = len(in_seqs)

    (tajimas_d, wattersons_theta) = tajimas_d_and_wattersons_theta(theta_pi, S, N)

    return((N, S, wattersons_theta, theta_pi, tajimas_d))


def calc_theta_pi_and_num_segregating(in_seqs):
    '''Calculate and return Theta Pi as well as the number of segregating
    sites as a tuple. Note: this function will not count polymorphisms
    where pairwise seqs do not have bases that are on of A, C, G, or T.'''

    bases = ['A', 'C', 'G', 'T']

    pairwise_comparisons = list(itertools.combinations(in_seqs, 2))

    segregating_sites = set()

    pairwise_diffs = []

    for combo in pairwise_comparisons:

        num_diff = 0

        seq1 = combo[0]
        seq2 = combo[1]

        if len(seq1) != len(seq2):
            print(seq1)
            print(seq2)
            sys.exit('Stopping - lengths of input sequences do not match. They should be aligned!')

        for i in range(len(seq1)):
        
            if seq1[i] != seq2[i]:
                if seq1[i] in bases and seq2[i] in bases:
                    num_diff += 1
                    segregating_sites.add(i)

        pairwise_diffs.append(num_diff)

    return((np.mean(pairwise_diffs), len(segregating_sites)))


def num_pairwise_diff_bases(bases):
    '''Will return the number of differences for a set of bases. Note: this
    function assumes that all bases are all one of A, C, T, G.'''

    num_diff = 0

    for possible_base in ['A', 'C', 'T', 'G']:
        
        num_specific_base = sum(x == possible_base for x in bases)
        num_other_base = len(bases) - num_specific_base

        if num_specific_base == 0 or num_other_base == 0:
            continue

        num_diff += num_specific_base * num_other_base
        bases = [b for b in bases if b != possible_base]

    return(num_diff)


def tajimas_d_and_wattersons_theta(theta_pi, num_seg_sites, n):
    '''Compute Tajima's D based on theta pi, number of segregating sites, and
    num_samples (n). Will return Watterson's theta as the second returned
    element as well.'''

    a1 = sum(1 / x for x in range(1, n))

    a2 = sum(1 / (x**2) for x in range(1, n))

    b1 = (n + 1) / (3 * (n - 1))

    b2 = (2 * (n**2 + n + 3)) / (9 * n * (n - 1))

    c1 = b1 - (1 / a1)

    c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / a1**2)  

    e1 = c1 / a1

    e2 = c2 / (a1**2 + a2)

    expected_sd = sqrt(e1 * num_seg_sites + e2 * num_seg_sites * (num_seg_sites - 1))

    wattersons_theta = num_seg_sites / a1

    if expected_sd == 0:
        tajimas_d = float("NaN")
    else:
        tajimas_d = (theta_pi - wattersons_theta) / expected_sd

    return(tajimas_d, wattersons_theta)

