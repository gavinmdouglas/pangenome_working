
rm(list = ls(all.names = TRUE))

library("discover")

species_presence <- read.table(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/summary/species_presence_core_25percent.tsv.gz",
                               header = TRUE, sep = "\t", row.names = 1)

Ellegaard.2019_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2019_SRR_ids.txt", stringsAsFactors = FALSE)$V1
Ellegaard.2020_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2020_SRR_ids.txt", stringsAsFactors = FALSE)$V1

species_presence <- t(species_presence)

species_presence <- species_presence[, c(Ellegaard.2019_samples, Ellegaard.2020_samples)]

species_strata_vec <- c(rep("Switzerland", length(Ellegaard.2019_samples)), rep("Japan", length(Ellegaard.2020_samples)))

species_presence <- species_presence[which(rowSums(species_presence) > 0), ]

species_events <- discover.matrix(species_presence, strata = species_strata_vec)

species_result.mutex <- pairwise.discover.test(species_events, fdr.method = "DBH", alternative = "less")
species_result.mutex_df <- as.data.frame(species_result.mutex, q.threshold = 1)

species_result.cooccur <- pairwise.discover.test(species_events, fdr.method = "DBH", alternative = "greater")
species_result.cooccur_df <- as.data.frame(species_result.cooccur, q.threshold = 1)

# Only significant hit is Gilliamella apicola vs. Bartonella apis.
length(which(as.numeric(species_presence["Gilliamella_apicola", ]) == 0 & as.numeric(species_presence["Bartonella_apis", ]) == 0))
length(which(as.numeric(species_presence["Gilliamella_apicola", ]) == 1 & as.numeric(species_presence["Bartonella_apis", ]) == 0))
length(which(as.numeric(species_presence["Gilliamella_apicola", ]) == 0 & as.numeric(species_presence["Bartonella_apis", ]) == 1))
length(which(as.numeric(species_presence["Gilliamella_apicola", ]) == 1 & as.numeric(species_presence["Bartonella_apis", ]) == 1))



# Strain presence / absence testing

species <- rownames(species_presence)

all_num_strains <- list()
for (sp in species) {
  AIC_file <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_gene_output/",
                    sp, "/", sp, ".strain_fit_summary.tsv", sep = "")
  AIC <- read.table(AIC_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  all_num_strains[[sp]] <- as.character(AIC[which(AIC$AIC == min(AIC$AIC)), "ID"])
}


all_strain_abun_raw <- list()

all_strains <- as.character()

for (sp in species) {
  sp_num_strains <- all_num_strains[[sp]]
  if (sp == "Bifidobacterium_asteroides") {
    strain_abun_file <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_gene_output/",
                              sp, "/best_fit_rerun/otu_table.", all_num_strains[[sp]], ".txt", sep = "")
  } else {
    strain_abun_file <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_gene_output/",
                              sp, "/strain_running/otu_table.", all_num_strains[[sp]], ".txt", sep = "")
  }
  strain_abun <- read.table(file = strain_abun_file, header = FALSE, sep = "\t", skip = 1)
  strain_abun[strain_abun < 0.01] <- 0
  strain_abun <- data.frame(sweep(strain_abun, 1, rowSums(strain_abun), '/')) * 100
  if (length(which(colSums(strain_abun) == 0)) > 0) {
    stop("Strain filtered out!") 
  }
  sp_sample_file <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_gene_input/prepped_input/",
                          sp, "_core_samples.txt", sep = "")
  sp_samples <- read.table(sp_sample_file, stringsAsFactors = FALSE)$V1
  rownames(strain_abun) <- sp_samples
  strain_abun[strain_abun > 0] <- 1
  
  colnames(strain_abun) <- paste(sp, as.character(1:ncol(strain_abun)), sep = ".")
  
  all_strain_abun_raw[[sp]] <- strain_abun
  
  all_strains <- c(all_strains, colnames(strain_abun))
}

strain_presence <- data.frame(matrix(0, nrow = length(all_strains), ncol = ncol(species_presence)))
colnames(strain_presence) <- colnames(species_presence)
rownames(strain_presence) <- all_strains

for (sp in species) {
  strain_abun <- data.frame(t(all_strain_abun_raw[[sp]]))
  strain_presence[rownames(strain_abun), colnames(strain_abun)] <- strain_abun
}

Ellegaard.2019_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2019_SRR_ids.txt", stringsAsFactors = FALSE)$V1
Ellegaard.2020_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2020_SRR_ids.txt", stringsAsFactors = FALSE)$V1

Ellegaard.2019_samples <- Ellegaard.2019_samples[which(Ellegaard.2019_samples %in% colnames(strain_presence))]
Ellegaard.2020_samples <- Ellegaard.2020_samples[which(Ellegaard.2020_samples %in% colnames(strain_presence))]

strain_presence <- strain_presence[, c(Ellegaard.2019_samples, Ellegaard.2020_samples)]

# Add species as rows of df as well.
species_strain_presence <- rbind(species_presence, strain_presence)

species_strain_strata_vec <- c(rep("Switzerland", length(Ellegaard.2019_samples)), rep("Japan", length(Ellegaard.2020_samples)))

species_strain_events <- discover.matrix(species_strain_presence, strata = species_strain_strata_vec)

species_strain_result.mutex <- pairwise.discover.test(species_strain_events, fdr.method = "DBH", alternative = "less")
species_strain_result.mutex_df <- as.data.frame(species_strain_result.mutex, q.threshold = 1)
species_strain_result.mutex_df$species1 <- gsub("\\..*$", "", species_strain_result.mutex_df$gene1)
species_strain_result.mutex_df$species2 <- gsub("\\..*$", "", species_strain_result.mutex_df$gene2)

species_strain_result.mutex_df_nonsamespecies <- species_strain_result.mutex_df[which(species_strain_result.mutex_df$species1 != species_strain_result.mutex_df$species2), ]
min(p.adjust(species_strain_result.mutex_df_nonsamespecies$p.value, "BH"))

species_strain_result.mutex_df_sig <- species_strain_result.mutex_df[which(species_strain_result.mutex_df$p.value < 0.05), ]
species_strain_result.mutex_df_sig <- species_strain_result.mutex_df_sig[which(species_strain_result.mutex_df_sig$species1 != species_strain_result.mutex_df_sig$species2), ]


species_strain_result.cooccur <- pairwise.discover.test(species_strain_events, fdr.method = "DBH", alternative = "greater")
species_strain_result.cooccur_df <- as.data.frame(species_strain_result.cooccur, q.threshold = 1)
species_strain_result.cooccur_df$species1 <- gsub("\\..*$", "", species_strain_result.cooccur_df$gene1)
species_strain_result.cooccur_df$species2 <- gsub("\\..*$", "", species_strain_result.cooccur_df$gene2)

species_strain_result.cooccur_df_nonsamespecies <- species_strain_result.cooccur_df[which(species_strain_result.cooccur_df$species1 != species_strain_result.cooccur_df$species2), ]
min(p.adjust(species_strain_result.cooccur_df_nonsamespecies$p.value, "BH"))


species_strain_result.cooccur_df_sig <- species_strain_result.cooccur_df[which(species_strain_result.cooccur_df$q.value < 0.2), ]
species_strain_result.cooccur_df_sig <- species_strain_result.cooccur_df_sig[which(species_strain_result.cooccur_df_sig$species1 != species_strain_result.cooccur_df_sig$species2), ]

# No strains were significantly mutually exclusive.
# Only strains within the same species were significantly co-occurring, which is trivial given how species are not found in all samples.
# Also, no strains of other species significantly co-occurring or co-exclusive with the presence of an overall species.
