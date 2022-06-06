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

sp_tajimasD_overall_raw <- list()
sp_tajimasD_samples_raw <- list()

tajimasD_working_dir <- "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/gene_haplotypes_tajimasD/"

for (sp in species_present) {
 
  tajimasD_summary_files <- list.files(path = paste(tajimasD_working_dir, sp, sep = "/"), full.names = TRUE, pattern = "_tajimas_d_and_metrics.tsv$") 
  
  sp_tajimasD_overall <- data.frame(matrix(NA, nrow = length(tajimasD_summary_files), ncol = 5))
  colnames(sp_tajimasD_overall) <- c("Species", "Gene", "Mean_samples", "Mean_all", "Mean_ref")
  sp_genes <- basename(gsub("_tajimas_d_and_metrics.tsv$", "", tajimasD_summary_files))
  rownames(sp_tajimasD_overall) <- sp_genes
  
  sp_tajimasD_samples <- data.frame(matrix(NA, nrow = length(tajimasD_summary_files), ncol = length(all_samples) + 2))
  rownames(sp_tajimasD_samples) <- sp_genes
  colnames(sp_tajimasD_samples) <- c("Species", "Gene", all_samples)
  
  sp_tajimasD_overall$Species <- sp
  sp_tajimasD_samples$Species <- sp
  sp_tajimasD_overall$Gene <- sp_genes
  sp_tajimasD_samples$Gene <- sp_genes
  
  for (i in 1:length(tajimasD_summary_files)) {
    raw_mean_tajimasD <- read.table(tajimasD_summary_files[i], header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
    
    tajimasD_samples_present <- rownames(raw_mean_tajimasD)[which(rownames(raw_mean_tajimasD) %in% all_samples)]
    
    sp_tajimasD_overall[i, c("Mean_samples", "Mean_all", "Mean_ref")] <- c(mean(raw_mean_tajimasD[tajimasD_samples_present, "D"]),
                                                                       raw_mean_tajimasD["all_haplotypes", "D"],
                                                                       raw_mean_tajimasD["reference_seqs", "D"])
    
    sp_tajimasD_samples[i, tajimasD_samples_present] <- raw_mean_tajimasD[tajimasD_samples_present, "D"]
  }
  
  sp_tajimasD_overall_raw[[sp]] <- sp_tajimasD_overall
  sp_tajimasD_samples_raw[[sp]] <- sp_tajimasD_samples

}

tajimasD_overall <- do.call(rbind, sp_tajimasD_overall_raw)
tajimasD_per_sample <- do.call(rbind, sp_tajimasD_samples_raw)

tajimasD_per_sample_long <- melt(tajimasD_per_sample, ids = c("Species", "Gene"),
                             variable.name = "Sample", value.name = "tajimasD")
tajimasD_per_sample_long <- tajimasD_per_sample_long[which(! is.na(tajimasD_per_sample_long$tajimasD)), ]


tajimasD_overall <- tajimasD_overall[which(tajimasD_overall$Mean_samples > 0), ]
tajimasD_overall$Species <- factor(tajimasD_overall$Species)


tajimasD_overall_long <- melt(tajimasD_overall, ids = c("Species", "Gene"), variable.name = "Type", value.name = "tajimasD")

tajimasD_overall$Mean_samples_vs_all <- tajimasD_overall$Mean_samples / tajimasD_overall$Mean_all

# Add in COG categories to each table.

tajimasD_per_sample_long$COG_category <- "No COG"
tajimasD_overall_long$COG_category <- "No COG"

tajimasD_per_sample_long_COG_genes_i <- which(tajimasD_per_sample_long$Gene %in% rownames(gene_COG_categories))
tajimasD_per_sample_long_COG_genes <- tajimasD_per_sample_long$Gene[tajimasD_per_sample_long_COG_genes_i]

tajimasD_overall_long_COG_genes_i <- which(tajimasD_overall_long$Gene %in% rownames(gene_COG_categories))
tajimasD_overall_long_COG_genes <- tajimasD_overall_long$Gene[tajimasD_overall_long_COG_genes_i]

tajimasD_per_sample_long$COG_category[tajimasD_per_sample_long_COG_genes_i] <- gene_COG_categories[tajimasD_per_sample_long_COG_genes, "COG_category"]
tajimasD_overall_long$COG_category[tajimasD_overall_long_COG_genes_i] <- gene_COG_categories[tajimasD_overall_long_COG_genes, "COG_category"]


tajimasD_overall_long <- tajimasD_overall_long[which(! is.na(tajimasD_overall_long$tajimasD)), ]
# Overall distributions by Species (based on mean Tajima's D values per gene)
ggplot(data = tajimasD_overall_long, aes(x = tajimasD, y = Species)) +
                                                geom_quasirandom(groupOnX=FALSE) +
                                                geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey") +
                                                scale_y_discrete(limits = rev(levels(tajimasD_overall$Species))) +
                                                theme_bw() +
                                                xlab("Tajima's D") +
                                                facet_grid(Type ~ .)

ggplot(data = tajimasD_overall, aes(x = Mean_samples_vs_all, y = Species)) +
                                            geom_quasirandom(groupOnX=FALSE) +
                                            geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey") +
                                            scale_x_continuous(trans='log10') +
                                            scale_y_discrete(limits = rev(levels(tajimasD_overall$Species))) +
                                            theme_bw() +
                                            xlab("Relative Tajima's D") +
                                            ggtitle("Ratio of Tajima's D based on alleles co-occurring in samples vs that of all alleles") +
                                            geom_vline(xintercept = 1, linetype = "dashed", colour = "yellow")




# Then look at how it varies by sample.
tajimasD_per_sample_long_simple <- tajimasD_per_sample_long[, -which(colnames(tajimasD_per_sample_long) %in% c("Gene", "Num_alleles", "COG_category"))]

# Remove species / sample combinations where there are fewer than 100 non-NA values.
for (species in unique(tajimasD_per_sample_long_simple$Species)) {
 
  for (samp in unique(tajimasD_per_sample_long_simple$Sample)) {
   
    species_samp_tajimasD <- tajimasD_per_sample_long_simple[which(tajimasD_per_sample_long_simple$Species == species & tajimasD_per_sample_long_simple$Sample == samp), "tajimasD"] 
    
    if (length(species_samp_tajimasD) == 0) { next }
    
    if (length(which(! is.na(species_samp_tajimasD))) < 100) {
      tajimasD_per_sample_long_simple <- tajimasD_per_sample_long_simple[-which(tajimasD_per_sample_long_simple$Species == species & tajimasD_per_sample_long_simple$Sample == samp), ]
    }
    
  }
   
}

# Set max to be 2 so that high tajimasD values can't skew distribution too much.
tajimasD_per_sample_long_simple$tajimasD[which(tajimasD_per_sample_long_simple$tajimasD > 2)] <- 2

tajimasD_per_sample_long_means <- aggregate(. ~ Species + Sample, data = tajimasD_per_sample_long_simple, FUN = function(x) { mean(x, na.rm = TRUE) })
tajimasD_per_sample_long_means$Species <- factor(tajimasD_per_sample_long_means$Species)

ggplot(data = tajimasD_per_sample_long_means, aes(x = tajimasD, y = Species)) +
  geom_quasirandom(groupOnX=FALSE) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey") +
  theme_bw() +
  xlab("Tajima's D") +
  ggtitle("Mean Tajima's D per sample and species (for instances with at least 100 datapoints)") +
  scale_y_discrete(limits = rev(levels(tajimasD_per_sample_long_means$Species)))


tajimasD_per_sample_long_medians <- aggregate(. ~ Species + Sample, data = tajimasD_per_sample_long_simple, FUN = function(x) { median(x, na.rm = TRUE) })
tajimasD_per_sample_long_medians$Species <- factor(tajimasD_per_sample_long_medians$Species)

ggplot(data = tajimasD_per_sample_long_medians, aes(x = tajimasD, y = Species)) +
  geom_quasirandom(groupOnX=FALSE) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey") +
  theme_bw() +
  xlab("Tajima's D") +
  ggtitle("Median Tajima's D per sample and species (for instances with at least 100 datapoints)") +
  scale_y_discrete(limits = rev(levels(tajimasD_per_sample_long_medians$Species)))


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


tajimasD_overall_long_multi <- expand_for_multi_COG_category(tajimasD_overall_long)

tajimasD_overall_long_multi_mean_samples <- tajimasD_overall_long_multi[which(tajimasD_overall_long_multi$Type == "Mean_samples"), ]


# Remove species / COG category combinations where there are fewer than 10 non-NA values.
for (species in unique(tajimasD_overall_long_multi_mean_samples$Species)) {
  
  for (category in unique(tajimasD_overall_long_multi_mean_samples$COG_category)) {
    
        species_COG_tajimasD <- tajimasD_overall_long_multi_mean_samples[which(tajimasD_overall_long_multi_mean_samples$Species == species & tajimasD_overall_long_multi_mean_samples$COG_category == category), "tajimasD"]
    
    if (length(species_COG_tajimasD) == 0) { next }
    
    if (length(which(! is.na(species_COG_tajimasD))) < 10) {
      tajimasD_overall_long_multi_mean_samples <- tajimasD_overall_long_multi_mean_samples[-which(tajimasD_overall_long_multi_mean_samples$Species == species & tajimasD_overall_long_multi_mean_samples$COG_category == category), ]
    }
    
  }
  
}

tajimasD_overall_long_multi_mean_samples$COG_category[which(tajimasD_overall_long_multi_mean_samples$COG_category == "No COG")] <- "-"


ggplot(data = tajimasD_overall_long_multi_mean_samples, aes(x = COG_category, y = tajimasD)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey") +
  theme_bw() +
  ylab("Tajima's D") +
  scale_y_continuous(trans='log10') +
  facet_wrap(Species ~ ., scales = "free_y")



# Look at how tajimasD varies by # alleles, freq in ref. genomes, freq. in MGS samples, and RVadj values
# Also for Bifidobacterium asteroides, look at Tajima's D by sig. non-sig.

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

tajimasD_overall_long_mean_samples <- tajimasD_overall_long[which(tajimasD_overall_long$Type == "Mean_samples"), ]

tajimasD_overall_long_mean_samples$Num_alleles <- num_alleles[tajimasD_overall_long_mean_samples$Gene, "Num_alleles"]

Num_alleles_lm_summary <- species_lm_summary_tab(in_df = tajimasD_overall_long_mean_samples, xvar = "Num_alleles", yvar = "log10(tajimasD + 0.00001)")

# Look at Tajima's D by species compared with the number of alleles
ggplot(data = tajimasD_overall_long_mean_samples, aes(x = Num_alleles, y = tajimasD)) +
  geom_point() +
  theme_bw() +
  ylab("Tajima's D (mean of samples)") +
  scale_y_continuous(trans='log10') +
  facet_wrap(Species ~ ., scales = "free_x") +
  geom_smooth(method = "lm")

mean(Num_alleles_lm_summary$beta, na.rm = TRUE)
mean(Num_alleles_lm_summary$adj_R_squared, na.rm = TRUE)


# Freq in ref. genomes
tajimasD_overall_long_mean_samples$ref_freq <- combined_panaroo[tajimasD_overall_long_mean_samples$Gene, "No..isolates"]

max_ref_freq <- c()
unique_sp <- unique(tajimasD_overall_long_mean_samples$Species)
for (sp in unique_sp) {
  max_ref_freq <- c(max_ref_freq, max(tajimasD_overall_long_mean_samples[which(tajimasD_overall_long_mean_samples$Species == sp), "ref_freq"]))
}
names(max_ref_freq) <- unique_sp
freq_ref_sp <- names(max_ref_freq)[which(max_ref_freq >= 10)]

tajimasD_overall_long_mean_samples_frequent_ref <- tajimasD_overall_long_mean_samples[which(tajimasD_overall_long_mean_samples$Species %in% freq_ref_sp), ]

ref_freq_lm_summary <- species_lm_summary_tab(in_df = tajimasD_overall_long_mean_samples_frequent_ref, xvar = "ref_freq", yvar = "log10(tajimasD + 0.00001)")

# Look at Tajima's D by species compared with the number of alleles
ggplot(data = tajimasD_overall_long_mean_samples_frequent_ref[which(tajimasD_overall_long_mean_samples_frequent_ref$Type == "Mean_samples"), ], aes(x = ref_freq, y = tajimasD)) +
  geom_point() +
  theme_bw() +
  ylab("Tajima's D (mean of samples)") +
  scale_y_continuous(trans='log10') +
  facet_wrap(Species ~ ., scales = "free_x") +
  geom_smooth(method = "lm")

mean(ref_freq_lm_summary$beta, na.rm = TRUE)
mean(ref_freq_lm_summary$adj_R_squared, na.rm = TRUE)

# Prevalence in MGS samples

tajimasD_overall_long_mean_samples$MGS_prev <- gene_prevalence[tajimasD_overall_long_mean_samples$Gene]

MGS_prev_lm_summary <- species_lm_summary_tab(in_df = tajimasD_overall_long_mean_samples, xvar = "MGS_prev", yvar = "log10(tajimasD + 0.00001)")

# Look at Tajima's D by species compared with the number of alleles
ggplot(data = tajimasD_overall_long_mean_samples, aes(x = MGS_prev, y = tajimasD)) +
  geom_point() +
  theme_bw() +
  ylab("Tajima's D (mean of samples)") +
  scale_y_continuous(trans='log10') +
  facet_wrap(Species ~ ., scales = "free_x") +
  geom_smooth(method = "lm")

mean(MGS_prev_lm_summary$beta, na.rm = TRUE)
mean(MGS_prev_lm_summary$adj_R_squared, na.rm = TRUE)


# RVadj values
RVadj_out_tested <- RVadj_out[which(! is.na(RVadj_out$mean_permuted_RVadj)), ]

tajimasD_overall_long_mean_samples_w_RVadj <- tajimasD_overall_long_mean_samples[which(tajimasD_overall_long_mean_samples$Gene %in% RVadj_out_tested$gene), ]

tajimasD_overall_long_mean_samples_w_RVadj$RVadj <- RVadj_out[tajimasD_overall_long_mean_samples_w_RVadj$Gene, "RVadj"]

tajimasD_overall_long_mean_samples_w_RVadj$RVadj_sig <- RVadj_out[tajimasD_overall_long_mean_samples_w_RVadj$Gene, "permuted_RVadj_p"] < 0.05

RVadj_lm_summary <- species_lm_summary_tab(in_df = tajimasD_overall_long_mean_samples_w_RVadj, xvar = "RVadj", yvar = "log10(tajimasD + 0.00001)")

# Look at Tajima's D by species compared with the number of alleles
ggplot(data = tajimasD_overall_long_mean_samples_w_RVadj, aes(x = RVadj, y = tajimasD)) +
  geom_point(aes(colour = RVadj_sig)) +
  theme_bw() +
  ylab("Tajima's D (mean of samples)") +
  scale_y_continuous(trans='log10') +
  facet_wrap(Species ~ ., scales = "free_x") +
  geom_smooth(method = "lm")

mean(RVadj_lm_summary$beta, na.rm = TRUE)
mean(RVadj_lm_summary$adj_R_squared, na.rm = TRUE)



# Try mixed model

# ROUGH CODE - LINEAR MODELS PER SPECIES MIGHT ACTUALLY MAKE THE MOST SENSE HERE, AS WE'RE INTERESTED IN HOW THE SPECIES MAY DIFFER TOO
# library(lme4)
# 
# tajimasD_overall_long_mean_samples_w_RVadj <- tajimasD_overall_long_mean_samples[which(! is.na(tajimasD_overall_long_mean_samples$RVadj)), ]
# 
# ref_freq_mixed_model_out <- lmer(log(tajimasD) ~ ref_freq + (1 | Species), data = tajimasD_overall_long_mean_samples_w_RVadj)
# ref_freq_MGS_prev_mixed_model_out <- lmer(log(tajimasD) ~ ref_freq + MGS_prev + (1 | Species), data = tajimasD_overall_long_mean_samples_w_RVadj)
# ref_freq_MGS_prev_Num_alleles_mixed_model_out <- lmer(log(tajimasD) ~ ref_freq + MGS_prev + Num_alleles + (1 | Species), data = tajimasD_overall_long_mean_samples_w_RVadj)
# ref_freq_MGS_prev_Num_alleles_RVadj_mixed_model_out <- lmer(log(tajimasD) ~ ref_freq + MGS_prev + Num_alleles + RVadj + (1 | Species), data = tajimasD_overall_long_mean_samples_w_RVadj)
# 
# mixed_model_out <- lmer(log(tajimasD) ~ ref_freq + MGS_prev + Num_alleles + RVadj + (1 | Species), data = tajimasD_overall_long_mean_samples_w_RVadj)
# summary(mixed_model_out)
# 
# anova(ref_freq_mixed_model_out, ref_freq_MGS_prev_mixed_model_out)
# anova(ref_freq_MGS_prev_mixed_model_out, ref_freq_MGS_prev_Num_alleles_mixed_model_out)
# anova(ref_freq_MGS_prev_Num_alleles_mixed_model_out, ref_freq_MGS_prev_Num_alleles_RVadj_mixed_model_out)


species_lm_output <- list()

for (sp in unique_sp) {
  tajimasD_overall_long_mean_samples_sp <- tajimasD_overall_long_mean_samples[which(tajimasD_overall_long_mean_samples$Species == sp), ]
  model_out <- lm(log10(tajimasD) ~ ref_freq + MGS_prev + Num_alleles, data = tajimasD_overall_long_mean_samples_sp)
  species_lm_output[[sp]] <- summary(model_out)
}

species_lm_output_rsquared <- sort(sapply(species_lm_output, function(x) { x$adj.r.squared }), decreasing = TRUE)
mean(species_lm_output_rsquared)


