#!/usr/bin/python

from __future__ import print_function

import argparse
import os
import sys
import vcf
import numpy as np

try:
   import cPickle as pickle
except:
   import pickle

def main():

    parser = argparse.ArgumentParser(

    description="Parse pandora-outputted VCF to numpy three-dimensional array (likely for input to StrainFinder).",

    epilog='''Usage example:

    python pandora_vcf_to_numpy_array.py --input VCF --genes GENES_FILE --output ARRAY

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="VCF", type=str,
                        help="Path to input Pandora VCF.", required=True)

    parser.add_argument("-g", "--genes", metavar="GENES_FILE", type=str,
                        help="Path to file that specified which genes to parse out - one gene per line.",
                        required=True)

    parser.add_argument("-s", "--samples", metavar="SAMPLES_FILE", type=str,
                        help="Path to file that specified which samples should be included (should be all those where the gene was called as present) - one per line.",
                        required=True)

    parser.add_argument("-m", "--min_depth", metavar="INT", type=int,
                        help="Total depth (across all alleles) to include a site (for a given sample).",
                        required=False, default=10)

    parser.add_argument("-c", "--cutoff", metavar="FLOAT", type=float,
                        help="Proportion of sites that must be called as present in order to include a sample.",
                        required=False, default=0.5)

    parser.add_argument("--keep_structural", action="store_true", required=False,
                        help="Flag to indicate that structural mutations (with a max of 4 alleles) will be output as well.")

    parser.add_argument("-o", "--output_prefix", metavar="ARRAY", type=str, required=True,
                        help="Prefix for output files (which will be the array, "
                             "and files indicating the order (and subset) of "
                             "samples and positions in the array. Will also "
                             "create an infofile with some statistics and "
                             "counts. If structural mutations are output "
                             "as well then the mapping of A, C, G, and T to "
                             "the actual structural variants at those "
                             "positions will also be output in a separate "
                             "file.")

    args = parser.parse_args()

    genes = []
    with open(args.genes, 'r') as gene_fh:
        for gene_line in gene_fh:
            genes.append(gene_line.rstrip())

    samples = []
    with open(args.samples, 'r') as sample_fh:
        for sample_line in sample_fh:
            samples.append(sample_line.rstrip())

    # Dict to keep track of structural alleles to be written to a separate file
    # (if --keep_structural was specified).
    structural_alleles = dict()
    structural_allele_ref = dict()

    # Counts to keep track of:
    num_struct_ignored = 0
    num_struct_too_many_alleles = 0
    total_intersecting_sites = 0
    num_skipped_due_to_preceeding_variant_at_same_site = 0

    sample_arrays = dict()
    for sample in samples:
        sample_arrays[sample] = np.empty((0, 4))

    variant_positions = []

    vcf_reader = vcf.Reader(open(args.input, 'r'))

    for gene in genes:

        try:
            tmp = vcf_reader.fetch(gene + ".fa")
        except:
            continue

        for record in vcf_reader.fetch(gene + ".fa"):
            
            scaffold = record.CHROM.replace(".fa", "")

            total_intersecting_sites += 1

            position = scaffold + ";" + str(record.POS)

            # Make sure all alleles are treated as strings and are uppercase.
            raw_alleles = record.alleles
            alleles = []
            for raw_allele in raw_alleles:
                alleles.append(str(raw_allele).upper())

            structural = False
            for allele in alleles:
                if len(allele) > 1:
                    structural = True
                    break

            if structural and not args.keep_structural:
                num_struct_ignored += 1
                continue

            num_samples = len(samples)
            num_starting_alleles = len(alleles)

            site_summary = np.arange(num_samples * num_starting_alleles).reshape(num_samples, num_starting_alleles)

            for i in range(num_samples):

                call = record.genotype(samples[i])

                allele_depth = np.array(call.data.MEAN_FWD_COVG) + np.array(call.data.MEAN_REV_COVG)

                # Set all depths to 0 for this sample if the total depth is below the cut-off.
                if allele_depth.sum() < args.min_depth:
                    allele_depth = np.zeros(num_starting_alleles)

                site_summary[i, :] = allele_depth

            # Only keep alleles that are present in any of these samples after the
            # filtering step.
            present_alleles = np.nonzero(site_summary.sum(axis=0))[0]

            if len(present_alleles) == 0:
                continue

            site_summary = site_summary[:, present_alleles]

            alleles = np.array(alleles)[present_alleles].tolist()

            num_alleles = len(alleles)

            #if num_alleles == 1 and alleles[0] == record.REF:
            #    continue

            if structural:

                # Ignore structural mutations with more than 4 classes.
                if num_alleles > 4:
                    num_struct_too_many_alleles += 1
                    continue
                elif num_alleles < 4:
                    num_missing_col = 4 - num_alleles

                    for i in range(num_missing_col):
                        site_summary = np.insert(site_summary, site_summary.shape[1],
                                                 np.zeros(num_samples), axis=1)

            else:
                # Need to fill in empty columns for SNPs if needed as well, and
                # make sure that allele columns are sorted alphabetically 
                site_summary_sorted = np.zeros((num_samples, 4))

                for i in range(len(alleles)):
                    site_summary_sorted[:, ['A', 'C', 'G', 'T'].index(alleles[i])] = site_summary[:, i]

                site_summary = site_summary_sorted

            # Record info for this site, unless there is already a variant recorded
            # for this site (which can happen sometimes in the pandora output)

            if position in variant_positions:
                num_skipped_due_to_preceeding_variant_at_same_site += 1
                continue

            if structural:
                structural_alleles[position] = alleles
                structural_allele_ref[position] = record.REF

            variant_positions.append(position)

            for sample_i in range(num_samples):
                sample = samples[sample_i]
                sample_arrays[sample] = np.insert(sample_arrays[sample],
                                                  sample_arrays[sample].shape[0],
                                                  site_summary[sample_i, :],
                                                  axis = 0)

    all_sample_data = []
    for sample in samples:
        all_sample_data.append(sample_arrays[sample])

    full_array = np.stack(all_sample_data, axis = 0)

    if full_array.shape[1] == 0:
        sys.exit("Stopping - no variant sites present so no output will be "
                 "produced.")

    full_array_depth_per_site = np.sum(full_array, axis=2)

    nonzero_prop = np.divide(np.count_nonzero(full_array_depth_per_site, axis=1), float(len(variant_positions)))

    samples_to_keep = list(np.argwhere(nonzero_prop > args.cutoff).flatten())

    num_samples_removed = len(samples) - len(samples_to_keep)

    full_array = full_array[samples_to_keep, :, :]

    # Then need to filter out any sites that are now missing in all samples.
    full_array_depth_per_site = full_array_depth_per_site[samples_to_keep, :]
    site_nonzero_counts = np.count_nonzero(full_array_depth_per_site, axis=0)
    sites_to_keep = list(np.argwhere(site_nonzero_counts > 0).flatten())

    full_array = full_array[:, sites_to_keep, :]

    sample_outfile = args.output_prefix + "_samples.txt"
    site_outfile = args.output_prefix + "_sites.txt"
    info_outfile = args.output_prefix + "_info.txt"
    pickle_outfile = args.output_prefix + ".np.cPickle"

    with open(sample_outfile, 'w') as sample_fh:
        for sample_i in samples_to_keep:
            print(samples[sample_i], file = sample_fh)

    if args.keep_structural:

        final_structural_pos = []
        with open(site_outfile, 'w') as site_fh:
            for site_i in sites_to_keep:
                print("\t".join(variant_positions[site_i].split(';')),
                      file = site_fh)

                if variant_positions[site_i] in structural_alleles.keys():
                    final_structural_pos.append(variant_positions[site_i])

        structured_alleles_outfile = args.output_prefix + "_structured_alleles.txt"
        with open(structured_alleles_outfile, 'w') as structured_alleles_fh:
            for structured_allele_pos in final_structural_pos:
                structured_allele_pos_clean = "\t".join(structured_allele_pos.split(';'))
                print(structured_allele_pos_clean + "\t" + structural_allele_ref[structured_allele_pos] + "\t" + ",".join(structural_alleles[structured_allele_pos]),
                      file = structured_alleles_fh)
    else:
        with open(site_outfile, 'w') as site_fh:
            for site_i in sites_to_keep:
                print(variant_positions[site_i], file = site_fh)

    with open(info_outfile, 'w') as info_fh:
        print("num_samples_called_as_present" + "\t" + str(len(samples)), file = info_fh)
        print("num_samples_removed" + "\t" + str(num_samples_removed), file = info_fh)
        print("num_final_samples" + "\t" + str(len(samples_to_keep)), file = info_fh)
        print("num_intersecting_sites_in_vcf" + "\t" + str(total_intersecting_sites), file = info_fh)
        print("num_final_sites" + "\t" + str(len(sites_to_keep)), file = info_fh)
        print("num_structured_sites_ignored" + "\t" + str(num_struct_ignored), file = info_fh)
        print("num_nonignored_structured_sites_w_too_many_alleles" + "\t" + str(num_struct_too_many_alleles), file = info_fh)
        print("num_skipped_due_to_preceeding_variant_at_same_site" + "\t" + str(num_skipped_due_to_preceeding_variant_at_same_site), file = info_fh)
        

    pickle.dump(full_array, open(pickle_outfile, 'w'))

if __name__ == '__main__':
    main()
