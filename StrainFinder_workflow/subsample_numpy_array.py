#!/usr/bin/python

from __future__ import print_function

import argparse
import os
import sys
import numpy as np
import random

try:
   import cPickle as pickle
except:
   import pickle

def main():

    parser = argparse.ArgumentParser(

        description="Subsample StrainFinder input 3-D array to the specified number of sites.",

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--input_array", metavar="ARRAY", type=str,
                        help="Path to input numpy array.", required=True)

    parser.add_argument("--input_sites", metavar="SITES", type=str,
                        help="Path to file that specified which genes to parse out - one gene per line.",
                        required=True)

    parser.add_argument("--num2keep", metavar="INT", type=int,
                        help="Number of sites to randomly subsample.",
                        required=True)

    parser.add_argument("-o", "--output_prefix", metavar="OUTPUT", type=str, required=True,
                        help="Prefix for output files (new array and sites files")

    args = parser.parse_args()

    inalign = pickle.load(open(args.input_array, 'rb'))

    num_sites = inalign.shape[1]

    if num_sites < args.num2keep:
        sys.exit("Error - fewer sites present than the number to subsample.")

    sites2keep = random.sample(range(num_sites), args.num2keep)

    inalign_subsampled = inalign[:, sites2keep, :]


    pickle_outfile = args.output_prefix + ".np.subsampled.cPickle"
    pickle.dump(inalign_subsampled, open(pickle_outfile, 'w'))

    sites_outfile = args.output_prefix + "_sites_subsampled.txt"
    sites_outfh = open(sites_outfile, 'w')

    line_i = 0
    with open(args.input_sites, 'r') as sites_infile:
        for line in sites_infile:
            if line_i in sites2keep:
                print(line, end = "", file = sites_outfh)
            line_i += 1

    sites_outfh.close()

if __name__ == '__main__':
    main()
