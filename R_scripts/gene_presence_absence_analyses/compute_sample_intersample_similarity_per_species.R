rm(list = ls(all.names = TRUE))

library("CooccurrenceAffinity")

# Compute inter-sample similarities based on non-cores genes per species

StrainFinder_core_gene_info <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/StrainFinder_core_genes_and_samples_inputs.rds")

gene_presence_noncore <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.noncore.matrix", header = TRUE, sep = "\t", row.names = 1)

species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species_present.txt", stringsAsFactors = FALSE, header = FALSE)$V1

per_sp_gene_presence_noncore <- list()
per_sp_similarities <- list()

for (sp in species) {
  
  sp_gene_presence_noncore <- gene_presence_noncore[grep(sp, rownames(gene_presence_noncore)), StrainFinder_core_gene_info[[sp]]$samples]
  
  sp_gene_presence_noncore <- sp_gene_presence_noncore[which(rowSums(sp_gene_presence_noncore) > 0), ]
  
  per_sp_gene_presence_noncore[[sp]] <- sp_gene_presence_noncore
  
  per_sp_similarities[[sp]] <- affinity(data = sp_gene_presence_noncore,
                                        row.or.col = "column",
                                        squarematrix = c("alpha_mle", "jaccard", "simpson"))
  
}

saveRDS(object = per_sp_gene_presence_noncore,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/gene_presence_noncore_by_species.rds")

saveRDS(object = per_sp_similarities,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/gene_presence_noncore_by_species_sample_similarities.rds")
