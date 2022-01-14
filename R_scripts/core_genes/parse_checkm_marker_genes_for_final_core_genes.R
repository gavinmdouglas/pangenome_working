### For species with < 10 genomes I identified core genes based on marker gene sets associated with the closest rank associated
### with the lineage in checkm. Need to figure out what panaroo orthologs these correspond to and whether they are found in the MGS data
### (at least in cases when most other potential core genes are also present).

rm(list = ls(all.names = TRUE))

species_subset <- c("Apilactobacillus_apinorum", "Apilactobacillus_kunkeei", "Bartonella_apis", 
                    "Bifidobacterium_coryneforme_indicum", "Bombella_sp", "Bombilactobacillus_mellifer",
                    "Frischella_perrara", "Gilliamella_sp", "Lactobacillus_apis",
                    "Lactobacillus_helsingborgensis", "Lactobacillus_kimbladii", "Lactobacillus_kullabergensis",
                    "Lactobacillus_melliventris", "Serratia_marcescens")

panaroo_passed_core <- list()

for (sp in species_subset) {

  panaroo_passed_core[[sp]] <- c()
  
  sp_checkm_marker_hits_file <- paste("/data1/gdouglas/projects/honey_bee/ref_genomes/checkm_core_genes/", sp, "_marker_hits.tsv", sep = "")
  sp_checkm_marker_hit_info <- read.table(sp_checkm_marker_hits_file, stringsAsFactors = FALSE, header = TRUE, sep = "\t")
  
  # Remove marker genes that have multiple hits in the same genome.
  sp_checkm_marker_hits <- sp_checkm_marker_hit_info$Gene.Id[grep(",", sp_checkm_marker_hit_info$Gene.Id, invert = TRUE)]
  
  sp_panaroo_outfile <- paste("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/", sp, "/gene_presence_absence_roary.csv.gz", sep = "")
  sp_panaroo <- read.table(sp_panaroo_outfile, header = TRUE, sep = ",", comment.char = "", quote = "")
  rownames(sp_panaroo) <- paste(sp, sp_panaroo$Gene, sep = "_")

  sp_num_genomes <- ncol(sp_panaroo) - 15 + 1
  
  colnames(sp_panaroo[, 15:ncol(sp_panaroo)])
  
  sp_potential_core <- list()
  
  for (gene in rownames(sp_panaroo)) {
    cds_ids <- sp_panaroo[gene, 15:ncol(sp_panaroo)]
    
    match_count = 0
    for (cds in cds_ids) {
       if (cds %in% sp_checkm_marker_hits) {
           match_count <- match_count + 1
       }
    }
    
    if (match_count == sp_num_genomes) {
      panaroo_passed_core[[sp]] <- c(panaroo_passed_core[[sp]], gene) 
    } else if(match_count > 0) {
      ### For sanity checks - choose some genes to check by hand.
      #print(gene) 
    }
    
  }
  
}


saveRDS(object = panaroo_passed_core,
        file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/draft_core_genes/checkm_panaroo_passed_potential_core.rds")

