#!/usr/bin/python3

from functions.iupac import check_ambig_match_dict

codon2aa = {"TTT":"F", "TTC":"F",
            "TTA":"L", "TTG":"L",
            "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
            "TAT":"Y", "TAC":"Y",
            "TAA":"STOP", "TAG":"STOP",
            "TGT":"C", "TGC":"C",
            "TGA":"STOP",
            "TGG":"W",
            "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
            "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
            "CAT":"H", "CAC":"H",
            "CAA":"Q", "CAG":"Q",
            "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
            "ATT":"I", "ATC":"I", "ATA":"I",
            "ATG":"M",
            "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
            "AAT":"N", "AAC":"N",
            "AAA":"K", "AAG":"K",
            "AGT":"S", "AGC":"S",
            "AGA":"R", "AGG":"R",
            "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
            "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
            "GAT":"D", "GAC":"D",
            "GAA":"E", "GAG":"E",
            "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}


def start_codon_present(seq):
    '''Check that start codon is present in an input DNA sequence of gene.
    Will call it present if it's within the first four codons. If any
    ambiguous bases present then will only call it not present if none of
    the possible codons match. Also, return the position of the start codon if
    present.'''

    start_codon_position = None

    expected_start_codon = seq[:3]

    if not check_ambig_match_dict(expected_start_codon, codon2aa, 'M'):

        exp_start_neighbour1 = seq[3:6]
        exp_start_neighbour2 = seq[6:9]
        exp_start_neighbour3 = seq[9:12]

        if check_ambig_match_dict(exp_start_neighbour1, codon2aa, 'M'):
            start_codon_present = True
            start_codon_position = 3
        elif check_ambig_match_dict(exp_start_neighbour2, codon2aa, 'M'):
            start_codon_present = True
            start_codon_position = 6
        elif check_ambig_match_dict(exp_start_neighbour3, codon2aa, 'M'):
            start_codon_present = True
            start_codon_position = 9
        else:
            start_codon_present = False
    else:
        start_codon_present = True
        start_codon_position = 0

    return((start_codon_present, start_codon_position))




def stop_codon_premature_present(seq, start_codon_position):
    '''Check whether premature stop codon is present in an input DNA sequence of
    gene. Also check if the expected stop codon is missing. Will call stop codon
    as present if it's within the last four codons (and stop codons within that
    range are not considered premature). If there are ambiguous bases in a codon
    then premature stop codons must be matched by all possible ambiguous codons,
    while the expected codons need only be matched by one. Also, takes in position
    of start codon, so that if the start codon is missing and/or a few codons in
    that a stop codon isn't called before it.'''

    if not start_codon_position:
        start_codon_position = 9

    premature_stop_codon_position = None
    percent_truncated = None

    stop_codon_calls = set()

    codon_positions = []

    for i in range(stop = len(seq), step = 3):

        codon = seq[i:i + 3]

        codon_positions.append(i)

        if i <= start_codon_position:
            continue

        if check_ambig_match_dict(codon, codon2aa, 'STOP',
                                  all_possible_match = True):
            stop_codon_calls.add(i)

    expected_start_codon = seq[:3]

    if not check_ambig_match_dict(expected_start_codon, codon2aa, 'M'):

        exp_start_neighbour1 = seq[3:6]
        exp_start_neighbour2 = seq[6:9]
        exp_start_neighbour3 = seq[9:12]

        if check_ambig_match_dict(exp_start_neighbour1, codon2aa, 'M'):
            start_codon_present = True
            start_codon_position = 3
        elif check_ambig_match_dict(exp_start_neighbour2, codon2aa, 'M'):
            start_codon_present = True
            start_codon_position = 6
        elif check_ambig_match_dict(exp_start_neighbour3, codon2aa, 'M'):
            start_codon_present = True
            start_codon_position = 9
        else:
            start_codon_present = False
    else:
        start_codon_present = True
        start_codon_position = 0

    return((start_codon_present, start_codon_position))

