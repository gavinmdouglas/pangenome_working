### Parse out read depths of core genes for species across all samples (for instances where species were identified as present).

rm(list = ls(all.names = TRUE))

species_presence <- read.table(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/summary/species_presence_core_90percent.tsv.gz",
                               header = TRUE, sep = "\t", row.names = 1)


bedgraph_files <- list.files("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/mean_depth_per_site", full.names = TRUE)
all_samples <- basename(gsub(".mean.bedGraph.gz", "", bedgraph_files))

tmp_in1 <- read.table(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/mean_depth_per_site/SRR10810002.mean.bedGraph.gz",
                      header = FALSE, sep = "\t", stringsAsFactors = FALSE)
gene_mean_depths <- data.frame(matrix(NA, nrow = nrow(tmp_in1), ncol = length(all_samples)))
colnames(gene_mean_depths) <- all_samples
rownames(gene_mean_depths) <- tmp_in1$V1

for (bedgraph in bedgraph_files) {
  samp <- basename(gsub(".mean.bedGraph.gz", "", bedgraph))
  tmp_in <- read.table(file = bedgraph, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  gene_mean_depths[tmp_in$V1, samp] <- tmp_in$V4
 }


# Save file for future reference.
# write.table(x = gene_mean_depths, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/summary/gene_mean_depth.tsv",
#             col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")


# Restrict to core genes.
panaroo_only_potential_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/panaroo_only_potential_core.rds")
checkm_panaroo_potential_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/checkm_panaroo_passed_potential_core.rds")
core_genes <- c(panaroo_only_potential_core, checkm_panaroo_potential_core)
all_core_genes <- as.character(unlist(core_genes))

genes_passed <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/ref_annot/all_species_pangenome_reference.trimmed.bed",
                           header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

all_core_genes <- all_core_genes[which(all_core_genes %in% genes_passed)]

core_gene_mean_depths <- gene_mean_depths[all_core_genes, ]

species_mean_depth <- species_presence
species_mean_depth[! is.na(species_mean_depth)] <- 0

for (sp in colnames(species_mean_depth)) {
 
  for (samp in rownames(species_mean_depth)) {
   
    if (species_presence[samp, sp] == 0) { next }
    
    sp_core_gene_depth <- core_gene_mean_depths[grep(sp, rownames(core_gene_mean_depths)), samp]
    
    species_mean_depth[samp, sp] <- mean(sp_core_gene_depth[which(sp_core_gene_depth > 0)])
     
  }
   
}


species_mean_depth <- species_mean_depth[, -which(colnames(species_mean_depth) %in% c("Apilactobacillus_apinorum", "Bombella_apis", "Bombella_sp"))]

species_mean_depth_rel <- species_mean_depth
species_mean_depth_rel <- data.frame(sweep(x = species_mean_depth, MARGIN = 1, STATS = rowSums(species_mean_depth), FUN = '/')) * 100


write.table(x = species_mean_depth, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/summary/species_mean_depth.tsv",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

write.table(x = species_mean_depth_rel, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/summary/species_rel_abun.tsv",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")
