#!/usr/bin/python3

from functions.codon import codon2aa
from math import log


def pairwise_dnds(seq1, seq2):
    
    possible_sites = exp_N_S_sites(seq1)

    obs_subs = obs_N_S_subs(seq1, seq2)

    pn = obs_subs[0] / possible_sites[0]

    ps = obs_subs[1] / possible_sites[1]

    dn = sub_prop_to_rate(pn)

    ds = sub_prop_to_rate(ps)

    if ds > 0:
        dnds = dn / ds
    else:
        dnds = float('NaN')        

    return((dn, ds, dnds))


def obs_N_S_subs(seq1, seq2):
    '''Number of nonsynonymous and synonymous substitutions between two sequences (codon-aligned already).'''
    
    if len(seq1) != len(seq2):
        sys.exit("Stopping - lengths of two input sequences differ.")

    Nd = 0
    Sd = 0

    for i in range(0, len(seq1) - 2, 3):
        
        codon1 = seq1[i:i + 3]
        codon2 = seq2[i:i + 3]

        if len(codon1) != 3 or '-' in codon1:
            continue

        if len(codon2) != 3 or '-' in codon2:
            continue

        if codon1 != codon2:

            non_standard_base = False
            bases = ['A', 'C', 'G', 'T']
            for codon_i in range(3):
                if codon1[codon_i] not in bases or codon2[codon_i] not in bases:
                    non_standard_base = True
                    break

            # Skip codons if any non-standard bases.
            if non_standard_base:
                continue

            # Note that there can be multiple subs in same codon so need to check each sub individually.
            for codon_i in range(3):
            
                codon1_base = codon1[codon_i]
                codon2_base = codon2[codon_i]

                if codon1_base != codon2_base:

                    new_codon = codon1[:codon_i] + codon2_base + codon2[codon_i + 1:]

                    codon1_aa = codon2aa[codon1]

                    if codon2aa[new_codon] == codon1_aa:
                        Sd += 1
                    else:
                        Nd += 1

    return((Nd, Sd))


def exp_N_S_sites(in_seq):
    '''Expected nonsynonymous and synonymous sites per sequence. Ignore codons
    that have any base besides four standard bases.'''

    bases = ['A', 'C', 'G', 'T']

    N = 0
    S = 0

    for i in range(0, len(in_seq) - 2, 3):
        
        codon = in_seq[i:i + 3]

        if len(codon) != 3 or "-" in codon:
            continue

        non_standard_base = False
        for codon_i in range(3):
            if codon[codon_i] not in bases:
                non_standard_base = True
                break

        # Skip codon if any non-standard bases.
        if non_standard_base:
            continue

        current_aa = codon2aa[codon]

        for codon_i in range(3):
            
            codon_base = codon[codon_i]

            for base in bases:

                if base != codon_base:

                    new_codon = codon[:codon_i] + base + codon[codon_i + 1:]

                    if codon2aa[new_codon] == current_aa:
                        S += 1 / 3
                    else:
                        N += 1 / 3

    return((N, S))


def sub_prop_to_rate(p):
    return(-1 * (3/4) * log(1 - (4 * p) / 3))

