### Code to run DISCOVER analysis between all genes  ###
### Run co-exclusion and co-occurence separately and also run 2019 and 2020 data together (stratified) ### 

rm(list = ls(all.names = TRUE))

library("discover")

gene_presence <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/summary/gene_presence_0.5_breadth.rds")

# Ignore all Apilactobacillus apinorum, Bombella apis and Bombella sp genes
gene_presence <- gene_presence[-grep("Apilactobacillus_apinorum", rownames(gene_presence)), ]
gene_presence <- gene_presence[-grep("Bombella_apis", rownames(gene_presence)), ]
gene_presence <- gene_presence[-grep("Bombella_sp", rownames(gene_presence)), ]

panaroo_only_potential_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/panaroo_only_potential_core.rds")
checkm_panaroo_potential_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/checkm_panaroo_passed_potential_core.rds")
core_genes <- c(panaroo_only_potential_core, checkm_panaroo_potential_core)

all_core_genes <- as.character(unlist(core_genes))

# Restrict to non-core genes.
gene_presence <- gene_presence[which(! rownames(gene_presence) %in% all_core_genes), ]

# Remove non-core genes present in more than 80% or fewer than 20% of samples
gene_prevalence <- rowSums(gene_presence > 0) / ncol(gene_presence)
gene_presence <- gene_presence[which(gene_prevalence > 0.8 | gene_prevalence < 0.2), ]


Ellegaard.2019_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2019_SRR_ids.txt", stringsAsFactors = FALSE)$V1
Ellegaard.2020_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2020_SRR_ids.txt", stringsAsFactors = FALSE)$V1

gene_presence <- gene_presence[, c(Ellegaard.2019_samples, Ellegaard.2020_samples)]

all_species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species.txt", stringsAsFactors = FALSE)$V1
gene_presence_species <- rep(NA, nrow(gene_presence))
for (sp in all_species) {
  gene_presence_species[grep(sp, rownames(gene_presence))] <- sp
}

table(gene_presence_species)


strata_vec <- c(rep("Switzerland", length(Ellegaard.2019_samples)), rep("Japan", length(Ellegaard.2020_samples)))

events <- discover.matrix(gene_presence, strata = strata_vec)

result.mutex <- pairwise.discover.test(events, fdr.method = "DBH", alternative = "less")
result.mutex_df <- as.data.frame(result.mutex, q.threshold = 1)

result.cooccur <- pairwise.discover.test(events, fdr.method = "DBH", alternative = "greater")
result.cooccur_df <- as.data.frame(result.cooccur, q.threshold = 1)

result.cooccur_df_DBH0.25 <- result.cooccur_df[which(result.cooccur_df$q.value < 0.25), ]

result.cooccur_df_DBH0.25$species1 <- NA
result.cooccur_df_DBH0.25$species2 <- NA
for (sp in all_species) {
  result.cooccur_df_DBH0.25$species1[grep(sp, result.cooccur_df_DBH0.25$gene1)] <- sp
  result.cooccur_df_DBH0.25$species2[grep(sp, result.cooccur_df_DBH0.25$gene2)] <- sp
}

#result.cooccur_df_DBH0.25_crossspecies <- result.cooccur_df_DBH0.25[which(result.cooccur_df_DBH0.25$species1 != result.cooccur_df_DBH0.25$species2), ]
saveRDS(object = result.cooccur_df_DBH0.25,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_bowtie2_mapping_DISCOVER_cooccur_results_clustered_DBH0.25.rds")
