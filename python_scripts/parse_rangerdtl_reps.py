#!/usr/bin/python3

import argparse
import os
import sys
import re
import numpy as np

def main():

    parser = argparse.ArgumentParser(

            description='Parse Ranger DTL replicate outputs to get mean number of duplications, transfers, and losses. Expects files to be in format PREFIX + REPLICATE #, e.g. PREFIX1, PREFIX2, etc.',

            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input', metavar='PREFIX', type=str,
                        help='Path to prefix for all input files.',
                        required=True)

    parser.add_argument('-n', '--name', metavar='NAME', type=str,
                        help='Gene name to print for output',
                        required=False, default = "gene")

    parser.add_argument('--max_reps', metavar='MAX_REPS', type=int,
                        help='Path to prefix for all input files.',
                        required=False, default=100)

    args = parser.parse_args()

    duplications = []
    transfers = []
    losses = []

    for rep_i in range(1, args.max_reps + 1, 1):

        rep_file = args.input + str(rep_i)

        with open(rep_file, "r") as INPUT:

            for line in INPUT:
                line = line.rstrip()
                line_split = line.split()

                if len(line) > 35 and line[0:35] == 'The minimum reconciliation cost is:':

                    p = re.compile(r'\d+')
                    all_matches = p.findall(line)

                    duplications.append(int(all_matches[1]))
                    transfers.append(int(all_matches[2]))
                    losses.append(int(all_matches[3]))


    print(args.name + '\t' + str(np.mean(duplications)) + '\t' + str(np.mean(transfers)) + '\t' + str(np.mean(losses)))

if __name__ == '__main__':
    main()
