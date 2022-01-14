### Call core genes for each species and make a single table with all genes mapped to.
### For species with < 10 genomes, choose core genes based on other data, e.g. checkm lineage core genes.
### For these other species a separate R script is used.

rm(list = ls(all.names = TRUE))

species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species.txt", stringsAsFactors = FALSE)$V1

species_to_ignore <- c("Apilactobacillus_apinorum", "Apilactobacillus_kunkeei", "Bartonella_apis", 
                      "Bifidobacterium_coryneforme_indicum", "Bombella_sp", "Bombilactobacillus_mellifer",
                      "Frischella_perrara", "Gilliamella_sp", "Lactobacillus_apis",
                      "Lactobacillus_helsingborgensis", "Lactobacillus_kimbladii", "Lactobacillus_kullabergensis",
                      "Lactobacillus_melliventris", "Serratia_marcescens")

core_genes <- list()

for (sp in species) {
  
  if (sp %in% species_to_ignore) { next }
   
  sp_panaroo_outfile <- paste("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/", sp, "/gene_presence_absence_roary.csv.gz", sep = "")
  
  sp_panaroo <- read.table(sp_panaroo_outfile, header = TRUE, sep = ",", comment.char = "", quote = "")
  
  rownames(sp_panaroo) <- paste(sp, sp_panaroo$Gene, sep = "_")
  
  sp_num_genomes <- ncol(sp_panaroo) - 15 + 1
  
  core_genes[[sp]] <- rownames(sp_panaroo)[which(sp_panaroo$No..isolates == sp_num_genomes)]
  
}

saveRDS(object = core_genes,
        file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/draft_core_genes/panaroo_only_potential_core.rds")
