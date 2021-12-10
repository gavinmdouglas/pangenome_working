#!/usr/bin/python3

import argparse
import os
import sys
import pandas as pd

def main():

    parser = argparse.ArgumentParser(

            description="Parse output directory and ensure that the correct files were created. Write out a descriptive message if any of this is not the case. Otherwise write out the inferred number of strains and the filepaths of the key files (that then would be easy to copy elsewhere).",

    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-g", "--gene_name", metavar="GENE_NAME", type=str,
                        help="Gene name (which is assumed to be the output folder name and prefix to output files).",
                        required=True)

    parser.add_argument('--parentfolder', metavar="PARENTFOLDER", type=str,
                        help="Folder that contains output folder the specific gene specified",
                        required=True)

    parser.add_argument('--subfolder', metavar="SUBFOLDER", type=str,
                        help="Subfolder of output folder that contains actual StrainFinder output files",
                        required=False, default='strain_running')

    parser.add_argument("--max_strains", metavar="MAX_STRAINS", type=int,
                        help="Num strains fit from 1 to max_strains",
                        required=False, default = 30)

    args = parser.parse_args()

    # Check if main output folder exists
    expected_path = args.parentfolder + '/' + args.gene_name
    if not os.path.exists(expected_path):
        print(args.gene_name + '\t' + 'no_folder')
        sys.exit()

    # If folder exists then check that StrainFinder subfolder exists
    strainfinder_expected_path = expected_path + '/' + args.subfolder
    if not os.path.exists(strainfinder_expected_path):
        print(args.gene_name + '\t' + 'no_strainfinder_subfolder')
        sys.exit()

    # If that exists then make sure that all expected output files are present.
    all_strainfinder_outputs = os.listdir(strainfinder_expected_path)
    for i in range(1, args.max_strains + 1, 1):
        em_file = strainfinder_expected_path + '/' + 'em.' + args.gene_name + '.' + str(i) + '.cpickle'
        log_file = strainfinder_expected_path + '/' + 'log.' + args.gene_name + '.' + str(i) + '.txt'
        otu_file = strainfinder_expected_path + '/' + 'otu_table.' + args.gene_name + '.' + str(i) + '.txt'

        if not os.path.exists(em_file) or not os.path.exists(log_file) or not os.path.exists(otu_file):
            print(args.gene_name + '\t' + 'only_partial_strainfinder_output')
            sys.exit()


    # Check that AIC table exists.
    AIC_file = expected_path + '/' + args.gene_name + '.strain_fit_summary.tsv'
    if not os.path.exists(AIC_file):
        print(args.gene_name + '\t' + 'no_aic_file')
        sys.exit()

    # Check that it's not malformed.
    AIC_in = pd.read_csv(AIC_file, sep = '\t')
    if AIC_in.columns[0] != 'ID' or AIC_in.columns[1] != 'AIC' or AIC_in.shape[0] != args.max_strains:
        print(args.gene_name + '\t' + 'malformed_AIC_table')
        sys.exit()

    # Get strain with lowest AIC
    inferred_num_strain = AIC_in.ID[AIC_in.AIC.idxmin()]


    # Make sure OTU table not malformed
    abun_file = strainfinder_expected_path + '/' + 'otu_table.' + args.gene_name + '.' + str(inferred_num_strain) + '.txt'
    otu_in = pd.read_csv(abun_file, sep = '\t')

    if otu_in.shape[1] != inferred_num_strain:
        print(args.gene_name + '\t' + 'malformed_OTU_table')
        sys.exit()

    # Passed to this point so print out info!
    print(args.gene_name + '\t' + 'passed' + '\t' + str(inferred_num_strain) + '\t' + abun_file + '\t' + AIC_file)


if __name__ == '__main__':
    main()
