#!/usr/bin/python3

import argparse
import os
import sys
import pysam
import pandas as pd
import numpy as np
from collections import defaultdict

class bed_coor:
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop


def main():

    parser = argparse.ArgumentParser(

            description='Parse merged BCF to output depth info needed for StrainFacts. Assumes that each separate contig is a separate gene that should be input separately into StrainFacts.',

    epilog='''Usage example:

    python bcfs_to_strainfacts_input.py --input INPUT --bed BED --outdir OUTDIR

    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input', metavar='INPUT', type=str, required=True,
                        help='Path to input BCF of all samples.')

    parser.add_argument('-b', '--bed', metavar='INPUT', type=str, required=True,
                        help='Path to reference bed file with gene names and coordinates of interest.')

    parser.add_argument('-o', '--output', metavar='OUTPUT', type=str,
                        help='Path to directory for output files.', required=True)

    args = parser.parse_args()

    bcf_in = pysam.VariantFile(args.input)

    gene_info = dict()
    metagenotypes = dict()

    with open(args.bed, 'r') as bed_filehandle:

        for bed_line in bed_filehandle:

            bed_line_split = bed_line.split()

            in_gene = bed_line_split[0]

            if in_gene != 'Bifidobacterium_asteroides_lspA':
                continue

            gene_info[in_gene] = bed_coor(int(bed_line_split[1]),
                                          int(bed_line_split[2]))

            num_pos = gene_info[in_gene].stop - gene_info[in_gene].start


    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    else:
        print('Stopping - output folder already exists.', file = sys.stderr)

#   for gene in gene_ids:

    gene = 'Bifidobacterium_asteroides_lspA'



    # Get sample ids, assuming that they are the first field after splitting by '.'
    # (and removing and preceding path if present).
    raw_samples = list(bcf_in.header.samples)
    samples = []
    for raw_sample in raw_samples:
        sample_basename = os.path.basename(raw_sample)
        samples.append(sample_basename.split('.')[0])





    SRR10810002_i = samples.index('SRR10810002')
    SRR10810029_i = samples.index('SRR10810029')
    SRR10810033_i = samples.index('SRR10810033')
    SRR7287247_i = samples.index('SRR7287247')





    num_samples = len(samples)

    site_info = []

    sample_alt_depth = defaultdict(list)
    sample_ref_depth = defaultdict(list)
    sample_total_depth = defaultdict(list)

    for rec in bcf_in.fetch(contig = gene,
                            start = gene_info[gene].start,
                            stop = gene_info[gene].stop):

        rec_split = str(rec).split()

        # Skip line if it is not polymorphic, or if it is multi-allelic.
        alt_allele = rec_split[4]
        if alt_allele == '.' or ',' in alt_allele:
            continue

        ref_allele = rec_split[3]

        site_pos = rec_split[1]

        site_info.append('_'.join([site_pos, ref_allele, alt_allele]))

        sample_genotypes = rec_split[-num_samples:]

        print('_'.join([site_pos, ref_allele, alt_allele]))
        print(sample_genotypes[SRR10810002_i])
        print(sample_genotypes[SRR10810029_i])
        print(sample_genotypes[SRR10810033_i])
        print(sample_genotypes[SRR7287247_i])
        print('\n\n')

        format_fields = rec_split[len(rec_split) - num_samples - 1].split(':')
        AD_index = format_fields.index('AD')

        for idx, geno_info in enumerate(sample_genotypes):

            allele_depth_raw = geno_info.split(':')[AD_index]
            allele_depth_raw = allele_depth_raw.replace('.', '0')
            allele_depth = allele_depth_raw.split(',')

            ref_depth = int(allele_depth[0])

            if len(allele_depth) == 2:
                alt_depth = int(allele_depth[1])
            elif len(allele_depth) == 1:
                alt_depth = 0
            else:
                print(geno_info, file = sys.stderr)
                sys.exit('Problem with allele depth?')

            sample_ref_depth[idx].append(ref_depth)
            sample_alt_depth[idx].append(alt_depth)
            sample_total_depth[idx].append(ref_depth + alt_depth)
    
    total_depth_df = pd.DataFrame.from_dict(sample_total_depth)

    sample_nonzero_sites_prop = total_depth_df.astype(bool).sum(axis=0) / total_depth_df.shape[0]
    total_depth_df_filt = total_depth_df.loc[:, sample_nonzero_sites_prop > 0.9]
    
    site_nonzero_sites_prop = total_depth_df.astype(bool).sum(axis=0) / total_depth_df.shape[0]
    total_depth_df_filt = total_depth_df_filt.loc[site_nonzero_sites_prop > 0.9, :]

    sample_subset = [samples[i] for i in list(total_depth_df_filt.columns)]
    site_info_subset = [site_info[i] for i in list(total_depth_df_filt.index)]

    if len(sample_subset) == 0 or len(site_info_subset) == 0:
        print('Skipping gene ' + gene + ' as at least all samples or all sites were dropped.', file = sys.stderr)
        next

    ref_depth_df = pd.DataFrame.from_dict(sample_ref_depth)
    alt_depth_df = pd.DataFrame.from_dict(sample_alt_depth)

    sample_ref_depth_filt = ref_depth_df.iloc[total_depth_df_filt.index, total_depth_df_filt.columns]
    sample_alt_depth_filt = alt_depth_df.iloc[total_depth_df_filt.index, total_depth_df_filt.columns]


    sample_ref_depth_filt.reset_index(inplace = True, drop = True)
    sample_ref_depth_filt.columns = range(sample_ref_depth_filt.shape[1])

    sample_ref_depth_filt_melt = pd.melt(sample_ref_depth_filt, ignore_index = False,
                                         value_name = 'metagenotype', var_name = 'sample')
    sample_ref_depth_filt_melt['position'] = sample_ref_depth_filt_melt.index
    sample_ref_depth_filt_melt['allele'] = 'ref'
    sample_ref_depth_filt_melt.reset_index(inplace = True, drop = True)


    sample_alt_depth_filt.reset_index(inplace = True, drop = True)
    sample_alt_depth_filt.columns = range(sample_alt_depth_filt.shape[1])

    sample_alt_depth_filt_melt = pd.melt(sample_alt_depth_filt, ignore_index = False,
                                         value_name = 'metagenotype', var_name = 'sample')
    sample_alt_depth_filt_melt['position'] = sample_alt_depth_filt_melt.index
    sample_alt_depth_filt_melt['allele'] = 'alt'
    sample_alt_depth_filt_melt.reset_index(inplace = True, drop = True)

    combined_gene_data = pd.concat([sample_ref_depth_filt_melt, sample_alt_depth_filt_melt],
                                   ignore_index = True, sort = False)

    combined_gene_data = combined_gene_data[['sample', 'position', 'allele', 'metagenotype']]

    outfile = args.output + '/' + gene + '_metagenotype.tsv'
    combined_gene_data.to_csv(path_or_buf = outfile, sep = '\t', index = False)

    sample_outfile = args.output + '/' + gene + '_samples.tsv'
    pd.Series(sample_subset).to_csv(path_or_buf = sample_outfile, sep = '\t', index = True, header = False)

    site_outfile = args.output + '/' + gene + '_sites.tsv'
    pd.Series(site_info_subset).to_csv(path_or_buf = site_outfile, sep = '\t', index = True, header = False)
            
if __name__ == '__main__':

    main()
