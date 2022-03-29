# Code should be put in an R notebook

rm(list = ls(all.names = TRUE))

library(ggbeeswarm)
library(ggplot2)
library(reshape2)
library(cowplot)

species_present <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species_present.txt",
                              stringsAsFactors = FALSE)$V1

all_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2019_and_2020_SRR_ids.txt",
                          stringsAsFactors = FALSE)$V1

num_alleles <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/num_alleles_per_gene.tsv",
                          header = FALSE, sep = "\t", row.names = 1)
colnames(num_alleles) <- c("Num_alleles")

gene_COG_categories <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/all_microbiota_COG.category_by_ortholog.tsv",
                                  header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
colnames(gene_COG_categories) <- "COG_category"

COG_category_descrip <- read.table("/data1/gdouglas/db/COG_definitions/COG_category_descrip.tsv",
                                   header = FALSE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
colnames(COG_category_descrip) <- "COG_category"

COG_category_descrip["No COG", "COG_category"] <- "No COG"


combined_panaroo <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/combined_panaroo_annot.rds")

gene_presence_absence <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.matrix",
                                    row.names = 1, header = TRUE, sep = "\t")
strain_abun <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/strain_abun.rds")
gene_prevalence <- c()
for (sp in names(strain_abun)) {
  sp_genes <- grep(sp, rownames(gene_presence_absence), value = TRUE)
  sp_genes_prev <- rowSums(gene_presence_absence[sp_genes, rownames(strain_abun[[sp]])])
  names(sp_genes_prev) <- sp_genes
  gene_prevalence <- c(gene_prevalence, sp_genes_prev)
}


RVadj_out <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/strain_vs_gene_abun_RVadj.rds")
rownames(RVadj_out) <- RVadj_out$gene

sp_dnds_overall_raw <- list()
sp_dnds_samples_raw <- list()

dnds_working_dir <- "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/gene_haplotypes_dNdS/"

for (sp in species_present) {
 
  dnds_summary_files <- list.files(path = paste(dnds_working_dir, sp, sep = "/"), full.names = TRUE, pattern = "_dnds.tsv$") 
  dnds_summary_files <- dnds_summary_files[grep("pairwise_haplotype_dnds.tsv$", dnds_summary_files, invert = TRUE)]
  
  sp_dnds_overall <- data.frame(matrix(NA, nrow = length(dnds_summary_files), ncol = 5))
  colnames(sp_dnds_overall) <- c("Species", "Gene", "Mean_samples", "Mean_all", "Mean_ref")
  sp_genes <- basename(gsub("_dnds.tsv", "", dnds_summary_files))
  rownames(sp_dnds_overall) <- sp_genes
  
  sp_dnds_samples <- data.frame(matrix(NA, nrow = length(dnds_summary_files), ncol = length(all_samples) + 2))
  rownames(sp_dnds_samples) <- sp_genes
  colnames(sp_dnds_samples) <- c("Species", "Gene", all_samples)
  
  sp_dnds_overall$Species <- sp
  sp_dnds_samples$Species <- sp
  sp_dnds_overall$Gene <- sp_genes
  sp_dnds_samples$Gene <- sp_genes
  
  for (i in 1:length(dnds_summary_files)) {
    raw_mean_dnds <- read.table(dnds_summary_files[i], header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
    
    dnds_samples_present <- rownames(raw_mean_dnds)[which(rownames(raw_mean_dnds) %in% all_samples)]
    
    sp_dnds_overall[i, c("Mean_samples", "Mean_all", "Mean_ref")] <- c(mean(raw_mean_dnds[dnds_samples_present, "dnds"]),
                                                                       raw_mean_dnds["all_haplotypes", "dnds"],
                                                                       raw_mean_dnds["reference_seqs", "dnds"])
    
    sp_dnds_samples[i, dnds_samples_present] <- raw_mean_dnds[dnds_samples_present, "dnds"]
  }
  
  sp_dnds_overall_raw[[sp]] <- sp_dnds_overall
  sp_dnds_samples_raw[[sp]] <- sp_dnds_samples

}

dnds_overall <- do.call(rbind, sp_dnds_overall_raw)
dnds_per_sample <- do.call(rbind, sp_dnds_samples_raw)

dnds_per_sample_long <- melt(dnds_per_sample, ids = c("Species", "Gene"),
                             variable.name = "Sample", value.name = "dN_dS")
dnds_per_sample_long <- dnds_per_sample_long[which(! is.na(dnds_per_sample_long$dN_dS)), ]


dnds_overall <- dnds_overall[which(dnds_overall$Mean_samples > 0), ]
dnds_overall$Species <- factor(dnds_overall$Species)


dnds_overall_long <- melt(dnds_overall, ids = c("Species", "Gene"), variable.name = "Type", value.name = "dN_dS")

dnds_overall$Mean_samples_vs_all <- dnds_overall$Mean_samples / dnds_overall$Mean_all

# Add in COG categories to each table.

dnds_per_sample_long$COG_category <- "No COG"
dnds_overall_long$COG_category <- "No COG"

dnds_per_sample_long_COG_genes_i <- which(dnds_per_sample_long$Gene %in% rownames(gene_COG_categories))
dnds_per_sample_long_COG_genes <- dnds_per_sample_long$Gene[dnds_per_sample_long_COG_genes_i]

dnds_overall_long_COG_genes_i <- which(dnds_overall_long$Gene %in% rownames(gene_COG_categories))
dnds_overall_long_COG_genes <- dnds_overall_long$Gene[dnds_overall_long_COG_genes_i]

dnds_per_sample_long$COG_category[dnds_per_sample_long_COG_genes_i] <- gene_COG_categories[dnds_per_sample_long_COG_genes, "COG_category"]
dnds_overall_long$COG_category[dnds_overall_long_COG_genes_i] <- gene_COG_categories[dnds_overall_long_COG_genes, "COG_category"]


# Overall distributions by Species (based on mean dN/dS values per gene)
ggplot(data = dnds_overall_long, aes(x = dN_dS, y = Species)) +
                                                geom_quasirandom(groupOnX=FALSE) +
                                                geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey") +
                                                scale_x_continuous(trans='log10') +
                                                scale_y_discrete(limits = rev(levels(dnds_overall$Species))) +
                                                theme_bw() +
                                                xlab("dN/dS") +
                                                facet_grid(Type ~ .)

ggplot(data = dnds_overall, aes(x = Mean_samples_vs_all, y = Species)) +
                                            geom_quasirandom(groupOnX=FALSE) +
                                            geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey") +
                                            scale_x_continuous(trans='log10') +
                                            scale_y_discrete(limits = rev(levels(dnds_overall$Species))) +
                                            theme_bw() +
                                            xlab("Relative dN/dS") +
                                            ggtitle("Ratio of dN/dS based on alleles co-occurring in samples vs that of all alleles") +
                                            geom_vline(xintercept = 1, linetype = "dashed", colour = "yellow")




# Then look at how it varies by sample.
dnds_per_sample_long_simple <- dnds_per_sample_long[, -which(colnames(dnds_per_sample_long) %in% c("Gene", "Num_alleles", "COG_category"))]

# Remove species / sample combinations where there are fewer than 100 non-NA values.
for (species in unique(dnds_per_sample_long_simple$Species)) {
 
  for (samp in unique(dnds_per_sample_long_simple$Sample)) {
   
    species_samp_dNdS <- dnds_per_sample_long_simple[which(dnds_per_sample_long_simple$Species == species & dnds_per_sample_long_simple$Sample == samp), "dN_dS"] 
    
    if (length(species_samp_dNdS) == 0) { next }
    
    if (length(which(! is.na(species_samp_dNdS))) < 100) {
      dnds_per_sample_long_simple <- dnds_per_sample_long_simple[-which(dnds_per_sample_long_simple$Species == species & dnds_per_sample_long_simple$Sample == samp), ]
    }
    
  }
   
}

# Set max to be 2 so that high dNdS values can't skew distribution too much.
dnds_per_sample_long_simple$dN_dS[which(dnds_per_sample_long_simple$dN_dS > 2)] <- 2

dnds_per_sample_long_means <- aggregate(. ~ Species + Sample, data = dnds_per_sample_long_simple, FUN = function(x) { mean(x, na.rm = TRUE) })
dnds_per_sample_long_means$Species <- factor(dnds_per_sample_long_means$Species)

ggplot(data = dnds_per_sample_long_means, aes(x = dN_dS, y = Species)) +
  geom_quasirandom(groupOnX=FALSE) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey") +
  theme_bw() +
  xlab("dN/dS") +
  ggtitle("Mean dN/dS per sample and species (for instances with at least 100 datapoints)") +
  scale_y_discrete(limits = rev(levels(dnds_per_sample_long_means$Species)))


dnds_per_sample_long_medians <- aggregate(. ~ Species + Sample, data = dnds_per_sample_long_simple, FUN = function(x) { median(x, na.rm = TRUE) })
dnds_per_sample_long_medians$Species <- factor(dnds_per_sample_long_medians$Species)

ggplot(data = dnds_per_sample_long_medians, aes(x = dN_dS, y = Species)) +
  geom_quasirandom(groupOnX=FALSE) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey") +
  theme_bw() +
  xlab("dN/dS") +
  ggtitle("Median dN/dS per sample and species (for instances with at least 100 datapoints)") +
  scale_y_discrete(limits = rev(levels(dnds_per_sample_long_medians$Species)))


# Look at breakdown per COG category
# This will take some pre-processing first, to make sure that genes mapped to multiple categories are treated right.

# Function to duplicate redundant rows that differ only be denoting the different COG categories a gene belongs to.
expand_for_multi_COG_category <- function(in_tab) {

  multi_row_i <- grep(",", in_tab$COG_category)
  
  if (length(multi_row_i) == 0) { return(in_tab) }
  
  multi_rows <- in_tab[multi_row_i, , drop = FALSE]
  final_tab <- in_tab[-multi_row_i, , drop = FALSE]
  
  new_rows <- list()
  new_row_index <- 1
  for (row_i in 1:nrow(multi_rows)) {
    row_content <-  multi_rows[row_i, ]
    
    assigned_categories <- strsplit(row_content$COG_category, ",")[[1]]
    
    for (assigned_category in assigned_categories) {
      tmp_row <- row_content
      tmp_row$COG_category <- assigned_category
      new_rows[[new_row_index]] <- tmp_row
      new_row_index <- new_row_index + 1
    }
  }
  
  new_multi <- do.call(rbind, new_rows)
  
  final_tab <- rbind(final_tab, new_multi)
  
  return(final_tab)
  
}


dnds_overall_long_multi <- expand_for_multi_COG_category(dnds_overall_long)

dnds_overall_long_multi_mean_samples <- dnds_overall_long_multi[which(dnds_overall_long_multi$Type == "Mean_samples"), ]


# Remove species / COG category combinations where there are fewer than 10 non-NA values.
for (species in unique(dnds_overall_long_multi_mean_samples$Species)) {
  
  for (category in unique(dnds_overall_long_multi_mean_samples$COG_category)) {
    
        species_COG_dNdS <- dnds_overall_long_multi_mean_samples[which(dnds_overall_long_multi_mean_samples$Species == species & dnds_overall_long_multi_mean_samples$COG_category == category), "dN_dS"]
    
    if (length(species_COG_dNdS) == 0) { next }
    
    if (length(which(! is.na(species_COG_dNdS))) < 10) {
      dnds_overall_long_multi_mean_samples <- dnds_overall_long_multi_mean_samples[-which(dnds_overall_long_multi_mean_samples$Species == species & dnds_overall_long_multi_mean_samples$COG_category == category), ]
    }
    
  }
  
}

dnds_overall_long_multi_mean_samples$COG_category[which(dnds_overall_long_multi_mean_samples$COG_category == "No COG")] <- "-"


ggplot(data = dnds_overall_long_multi_mean_samples, aes(x = COG_category, y = dN_dS)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey") +
  theme_bw() +
  ylab("dN/dS") +
  scale_y_continuous(trans='log10') +
  facet_wrap(Species ~ ., scales = "free_y")



# Look at how dnds varies by # alleles, freq in ref. genomes, freq. in MGS samples, and RVadj values
# Also for Bifidobacterium asteroides, look at dN/dS by sig. non-sig.

species_lm_summary_tab <- function(in_df, xvar, yvar) {
  
  all_species <- unique(in_df$Species)
  
  out_df <- data.frame(matrix(NA, nrow = length(all_species), ncol = 4))
  rownames(out_df) <- all_species
  colnames(out_df) <- c("intercept", "beta", "adj_R_squared", "p")
  
  for (sp in all_species) {

    sp_lm_out <- lm(as.formula(paste(yvar, "~", xvar)),
                    data = in_df[which(in_df$Species == sp), ])
    sp_lm_out_summary <- summary(sp_lm_out)
  
    out_df[sp, "p"] <- sp_lm_out_summary$coefficients[2, 4]
    
    if (out_df[sp, "p"] < 0.05) {
    
      out_df[sp, c("intercept", "beta")] <- sp_lm_out$coefficients
      out_df[sp, "adj_R_squared"] <- sp_lm_out_summary$adj.r.squared
    
    }
    
  }
  
  return(out_df)
  
}


# Num_alleles

dnds_overall_long_mean_samples <- dnds_overall_long[which(dnds_overall_long$Type == "Mean_samples"), ]

dnds_overall_long_mean_samples$Num_alleles <- num_alleles[dnds_overall_long_mean_samples$Gene, "Num_alleles"]

Num_alleles_lm_summary <- species_lm_summary_tab(in_df = dnds_overall_long_mean_samples, xvar = "Num_alleles", yvar = "log10(dN_dS + 0.00001)")

# Look at dN/dS by species compared with the number of alleles
ggplot(data = dnds_overall_long_mean_samples, aes(x = Num_alleles, y = dN_dS)) +
  geom_point() +
  theme_bw() +
  ylab("dN/dS (mean of samples)") +
  scale_y_continuous(trans='log10') +
  facet_wrap(Species ~ ., scales = "free_x") +
  geom_smooth(method = "lm")

mean(Num_alleles_lm_summary$beta, na.rm = TRUE)
mean(Num_alleles_lm_summary$adj_R_squared, na.rm = TRUE)


# Freq in ref. genomes
dnds_overall_long_mean_samples$ref_freq <- combined_panaroo[dnds_overall_long_mean_samples$Gene, "No..isolates"]

max_ref_freq <- c()
unique_sp <- unique(dnds_overall_long_mean_samples$Species)
for (sp in unique_sp) {
  max_ref_freq <- c(max_ref_freq, max(dnds_overall_long_mean_samples[which(dnds_overall_long_mean_samples$Species == sp), "ref_freq"]))
}
names(max_ref_freq) <- unique_sp
freq_ref_sp <- names(max_ref_freq)[which(max_ref_freq >= 10)]

dnds_overall_long_mean_samples_frequent_ref <- dnds_overall_long_mean_samples[which(dnds_overall_long_mean_samples$Species %in% freq_ref_sp), ]

ref_freq_lm_summary <- species_lm_summary_tab(in_df = dnds_overall_long_mean_samples_frequent_ref, xvar = "ref_freq", yvar = "log10(dN_dS + 0.00001)")

# Look at dN/dS by species compared with the number of alleles
ggplot(data = dnds_overall_long_mean_samples_frequent_ref[which(dnds_overall_long_mean_samples_frequent_ref$Type == "Mean_samples"), ], aes(x = ref_freq, y = dN_dS)) +
  geom_point() +
  theme_bw() +
  ylab("dN/dS (mean of samples)") +
  scale_y_continuous(trans='log10') +
  facet_wrap(Species ~ ., scales = "free_x") +
  geom_smooth(method = "lm")

mean(ref_freq_lm_summary$beta, na.rm = TRUE)
mean(ref_freq_lm_summary$adj_R_squared, na.rm = TRUE)

# Prevalence in MGS samples

dnds_overall_long_mean_samples$MGS_prev <- gene_prevalence[dnds_overall_long_mean_samples$Gene]

MGS_prev_lm_summary <- species_lm_summary_tab(in_df = dnds_overall_long_mean_samples, xvar = "MGS_prev", yvar = "log10(dN_dS + 0.00001)")

# Look at dN/dS by species compared with the number of alleles
ggplot(data = dnds_overall_long_mean_samples, aes(x = MGS_prev, y = dN_dS)) +
  geom_point() +
  theme_bw() +
  ylab("dN/dS (mean of samples)") +
  scale_y_continuous(trans='log10') +
  facet_wrap(Species ~ ., scales = "free_x") +
  geom_smooth(method = "lm")

mean(MGS_prev_lm_summary$beta, na.rm = TRUE)
mean(MGS_prev_lm_summary$adj_R_squared, na.rm = TRUE)


# RVadj values
RVadj_out_tested <- RVadj_out[which(! is.na(RVadj_out$mean_permuted_RVadj)), ]

dnds_overall_long_mean_samples_w_RVadj <- dnds_overall_long_mean_samples[which(dnds_overall_long_mean_samples$Gene %in% RVadj_out_tested$gene), ]

dnds_overall_long_mean_samples_w_RVadj$RVadj <- RVadj_out[dnds_overall_long_mean_samples_w_RVadj$Gene, "RVadj"]

dnds_overall_long_mean_samples_w_RVadj$RVadj_sig <- RVadj_out[dnds_overall_long_mean_samples_w_RVadj$Gene, "permuted_RVadj_p"] < 0.05

RVadj_lm_summary <- species_lm_summary_tab(in_df = dnds_overall_long_mean_samples_w_RVadj, xvar = "RVadj", yvar = "log10(dN_dS + 0.00001)")

# Look at dN/dS by species compared with the number of alleles
ggplot(data = dnds_overall_long_mean_samples_w_RVadj, aes(x = RVadj, y = dN_dS)) +
  geom_point(aes(colour = RVadj_sig)) +
  theme_bw() +
  ylab("dN/dS (mean of samples)") +
  scale_y_continuous(trans='log10') +
  facet_wrap(Species ~ ., scales = "free_x") +
  geom_smooth(method = "lm")

mean(RVadj_lm_summary$beta, na.rm = TRUE)
mean(RVadj_lm_summary$adj_R_squared, na.rm = TRUE)



# Try mixed model

# ROUGH CODE - LINEAR MODELS PER SPECIES MIGHT ACTUALLY MAKE THE MOST SENSE HERE, AS WE'RE INTERESTED IN HOW THE SPECIES MAY DIFFER TOO
# library(lme4)
# 
# dnds_overall_long_mean_samples_w_RVadj <- dnds_overall_long_mean_samples[which(! is.na(dnds_overall_long_mean_samples$RVadj)), ]
# 
# ref_freq_mixed_model_out <- lmer(log(dN_dS) ~ ref_freq + (1 | Species), data = dnds_overall_long_mean_samples_w_RVadj)
# ref_freq_MGS_prev_mixed_model_out <- lmer(log(dN_dS) ~ ref_freq + MGS_prev + (1 | Species), data = dnds_overall_long_mean_samples_w_RVadj)
# ref_freq_MGS_prev_Num_alleles_mixed_model_out <- lmer(log(dN_dS) ~ ref_freq + MGS_prev + Num_alleles + (1 | Species), data = dnds_overall_long_mean_samples_w_RVadj)
# ref_freq_MGS_prev_Num_alleles_RVadj_mixed_model_out <- lmer(log(dN_dS) ~ ref_freq + MGS_prev + Num_alleles + RVadj + (1 | Species), data = dnds_overall_long_mean_samples_w_RVadj)
# 
# mixed_model_out <- lmer(log(dN_dS) ~ ref_freq + MGS_prev + Num_alleles + RVadj + (1 | Species), data = dnds_overall_long_mean_samples_w_RVadj)
# summary(mixed_model_out)
# 
# anova(ref_freq_mixed_model_out, ref_freq_MGS_prev_mixed_model_out)
# anova(ref_freq_MGS_prev_mixed_model_out, ref_freq_MGS_prev_Num_alleles_mixed_model_out)
# anova(ref_freq_MGS_prev_Num_alleles_mixed_model_out, ref_freq_MGS_prev_Num_alleles_RVadj_mixed_model_out)


species_lm_output <- list()

for (sp in unique_sp) {
  dnds_overall_long_mean_samples_sp <- dnds_overall_long_mean_samples[which(dnds_overall_long_mean_samples$Species == sp), ]
  model_out <- lm(log10(dN_dS) ~ ref_freq + MGS_prev + Num_alleles, data = dnds_overall_long_mean_samples_sp)
  species_lm_output[[sp]] <- summary(model_out)
}

species_lm_output_rsquared <- sort(sapply(species_lm_output, function(x) { x$adj.r.squared }), decreasing = TRUE)



