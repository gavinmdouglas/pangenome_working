### Explore breadth of coverage in breadth of coverage output.
### Call species as present across specific samples based on core gene presence profiles.

rm(list = ls(all.names = TRUE))

library(ggplot2)

all_present <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/summary/gene_presence_0.5_breadth.rds")
breadth_output <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/summary/gene_breadth_breakdown.rds")

panaroo_only_potential_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/panaroo_only_potential_core.rds")
checkm_panaroo_potential_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/checkm_panaroo_passed_potential_core.rds")
core_genes <- c(panaroo_only_potential_core, checkm_panaroo_potential_core)

genes_passed <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/ref_annot/all_species_pangenome_reference.trimmed.bed",
                                header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1
all_species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species.txt", stringsAsFactors = FALSE)$V1
# Remove core genes that did not make it above the size cut-off.
for (sp in all_species) {
  core_genes[[sp]] <- core_genes[[sp]][which(core_genes[[sp]] %in% genes_passed)]
}

# Distribution of core and non-core genes.
ggplot(data = breadth_output$breadth_summary, aes(y = species, x = num_at_least_0.5)) +
       geom_boxplot() +
       facet_wrap(gene_type ~ .) +
       scale_y_discrete(limits = rev) +
       xlab("No. samples where breadth >= 0.5") +
       ylab("Species")




# Breakdown of % core genes (of those above size cut-off) called as present per sample per species.

all_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2019_and_2020_SRR_ids.txt", stringsAsFactors = FALSE)$V1

core_gene_breakdown <- data.frame(matrix(NA, nrow = length(all_species) * length(all_samples), ncol = 5))
colnames(core_gene_breakdown) <- c("Species", "Sample", "Num_any", "Num_core", "Percent_core")

core_gene_breakdown$Species <- rep(x = all_species, each = length(all_samples))
core_gene_breakdown$Sample <- rep(x = all_samples, times = length(all_species))
  
for (sp in all_species) {
  
  sp_core <- core_genes[[sp]]
  
  sp_row_i <- grep(sp, rownames(all_present))
  
  for (samp in all_samples) {
 
    row_i <- which(core_gene_breakdown$Species == sp & core_gene_breakdown$Sample == samp)
    
    core_gene_breakdown[row_i, "Num_any"] <- length(which(all_present[sp_row_i, samp] > 0))
    
    core_gene_breakdown[row_i, "Num_core"] <- length(which(all_present[sp_core, samp] > 0))
    
    core_gene_breakdown[row_i, "Percent_core"] <- (core_gene_breakdown[row_i, "Num_core"] / length(sp_core)) * 100
    
  }

}


ggplot(data = core_gene_breakdown, aes(y = Species, x = Percent_core)) +
  geom_boxplot() +
  scale_y_discrete(limits = rev) +
  xlab("Percent core genes present") +
  ylab("Species") +
  xlim(0, 100)



# Scatterplot of no. core genes vs all species genes.
ggplot(data = core_gene_breakdown, aes(x = Num_any, y = Percent_core)) +
  geom_point() +
  facet_wrap(Species ~ .)
  


# Get species presence profile based on 25% and 90% core genes present.
species_25per <- data.frame(matrix(0, nrow = length(all_samples), ncol = length(all_species)))
colnames(species_25per) <- all_species
rownames(species_25per) <- all_samples

species_90per <- species_25per

for (sp in all_species) {
  
  sp_core <- core_genes[[sp]]
  
  for (samp in all_samples) {
 
    percent_core <- core_gene_breakdown[which(core_gene_breakdown$Species == sp & core_gene_breakdown$Sample == samp), "Percent_core"]
    
    if (percent_core >= 25) {
      species_25per[samp, sp] <- 1
    }
    
    if (percent_core >= 90) {
      species_90per[samp, sp] <- 1
    }
    
  }
  
}

write.table(x = species_25per, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/summary/species_presence_core_25percent.tsv",
            row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")

write.table(x = species_90per, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/summary/species_presence_core_90percent.tsv",
            row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
