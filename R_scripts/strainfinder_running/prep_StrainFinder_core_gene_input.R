# Write out sample and gene files for running StrainFinder on core genes to infer strains.

rm(list = ls(all.names = TRUE))

StrainFinder_core_gene_info <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/StrainFinder_core_genes_and_samples_inputs.rds")

for (sp in names(StrainFinder_core_gene_info)) {
  sample_outfile <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_gene_input/starting_files/",
                          sp, "_samples.txt", sep = "")
  
  gene_outfile <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_gene_input/starting_files/",
                          sp, "_genes.txt", sep = "")
  
  write.table(x = StrainFinder_core_gene_info[[sp]]$samples, file = sample_outfile, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(x = StrainFinder_core_gene_info[[sp]]$genes, file = gene_outfile, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}
