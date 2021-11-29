import cPickle
import os
import sys
import vcf
import numpy as np
import pandas as pd

all_genes_w_variation = []
all_genes_w_variation_filepath = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_individual_genes/all_genes_w_variation.txt"

with open(all_genes_w_variation_filepath, 'r') as gene_id_file:
    for gene in gene_id_file:
        all_genes_w_variation.append(gene.rstrip())

all_genes_w_variation = np.asarray(all_genes_w_variation)


all_samples = []
all_samples_filepath = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2019_and_2020_SRR_ids.txt"
with open(all_samples_filepath, 'r') as all_samples_file:
    for sample in all_samples_file:
        all_samples.append(sample.rstrip())
all_samples = np.asarray(all_samples)


mean_depth_of_variants = pd.DataFrame(np.empty((len(all_genes_w_variation), len(all_samples))),
                                      index = all_genes_w_variation,
                                      columns=all_samples)


mean_depth_of_variants.loc[:, :] = 0

for gene in all_genes_w_variation:

    samples_filepath = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_individual_genes/prepped_input/" + gene + "_samples.txt"
    samples = []
    with open(samples_filepath, 'r') as sample_file:
        for sample in sample_file:
            samples.append(sample.rstrip())

    pickle_filepath = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_individual_genes/prepped_input/" + gene + ".np.cPickle"
    inalign = cPickle.load(open(pickle_filepath, 'rb'))

    for s_i in range(len(samples)):

        depth_per_site = []

        for site_i in range(inalign.shape[1]):

            depth_per_site.append(inalign[s_i, site_i, 0] + inalign[s_i, site_i, 1] + inalign[s_i, site_i, 2] + inalign[s_i, site_i, 3])

        mean_depth_of_variants.loc[gene, samples[s_i]] = sum(depth_per_site) / len(depth_per_site)


mean_depth_of_variants.to_csv(path_or_buf = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina/pandora_mean_variant_site_depth.tsv",
                              sep = '\t')

