#!/usr/bin/python3

import argparse
import os
import numpy as np
import pandas as pd

def main():

    parser = argparse.ArgumentParser(

            description="Compute mean depth per sample per separate "
                        "chrom/scaffold in a samtools depth output table. "
                        "This assumes that the table has a header line "
                        "(output with the -H option).",

    epilog='''Usage example:

    python mean_samtools_depth.py --input TABLE --clean_colnames \
           --output OUTFILE

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="PATH", type=str,
                        help="Path to input table corresponding to output of "
                             "samtools depth (with the -H option used to get "
                             "a headerline)", required=True)

    parser.add_argument("-o", "--output", metavar="OUTFILE", type=str,
                        help="Path to output table", required=True)

    parser.add_argument("-c", "--clean_colnames", action='store_true',
                        help="Set to clean up column names to be just file "
                             "basenames (no path or suffix)",
                        required = False, default = False)

    args = parser.parse_args()


    depth_table = pd.read_table(filepath_or_buffer = args.input,
                                sep = "\t",
                                index_col = ["#CHROM", "POS"],
                                header = 0)

    new_columns = []
    for colname in depth_table.columns:
        new_columns.append(os.path.splitext(os.path.basename(colname))[0])

    if args.clean_colnames:
        depth_table.columns = new_columns

    depth_table = depth_table.groupby(level=0).mean()

    depth_table.to_csv(path_or_buf=args.output, sep='\t', header=True,
                       index=True)


if __name__ == '__main__':
    main()
