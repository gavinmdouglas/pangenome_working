rm(list = ls(all.names = TRUE))

library(ggbeeswarm)
library(ggplot2)
library(reshape2)
library(cowplot)

source("/home/gdouglas/scripts/pangenome_working/R_scripts/functions.R")

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

tajimasD_working_dir <- "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/gene_haplotypes_tajimasD/"

for (sp in species_present) {
 
  tajimasD_summary_files <- list.files(path = paste(tajimasD_working_dir, sp, sep = "/"), full.names = TRUE, pattern = "tajimas_d_and_metrics.tsv$") 

  sp_tajimasD_overall <- data.frame(matrix(NA, nrow = length(tajimasD_summary_files), ncol = 7))
  colnames(sp_tajimasD_overall) <- c("Species", "Gene", "mean_n", "mean_S", "mean_Wattersons", "mean_pi", "mean_D")
  sp_genes <- basename(gsub("_tajimas_d_and_metrics.tsv", "", tajimasD_summary_files))
  rownames(sp_tajimasD_overall) <- sp_genes
  
  sp_tajimasD_overall$Species <- sp
  sp_tajimasD_overall$Gene <- sp_genes
  
  for (i in 1:length(tajimasD_summary_files)) {
    raw_mean_tajimasD <- read.table(tajimasD_summary_files[i], header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
    sp_tajimasD_overall[i, c("mean_n", "mean_S", "mean_Wattersons", "mean_pi", "mean_D")] <- raw_mean_tajimasD["all_samples", ]
  }
  
  sp_tajimasD_overall_raw[[sp]] <- sp_tajimasD_overall

}

tajimasD_overall <- do.call(rbind, sp_tajimasD_overall_raw)

tajimasD_overall$Species <- factor(tajimasD_overall$Species)
tajimasD_overall$num_alleles <- num_alleles[tajimasD_overall$Gene, ]



# Exclude genes with fewer than 6 alleles (as a permutation test based on that many allelic combos is almost certainly not useful).
# Chosen because there are 720 different permutations of 6 elements, whereas there are only 120 permutations of 5 elements.
tajimasD_overall_norare <- tajimasD_overall[which(tajimasD_overall$num_alleles >= 6), ]


# Also exclude core genes, as these are especially linked to overall strain distances.
gene_presence_noncore <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.noncore.matrix", header = TRUE, sep = "\t", row.names = 1)

tajimasD_overall_norare_noncore <- tajimasD_overall_norare[which(tajimasD_overall_norare$Gene %in% rownames(gene_presence_noncore)), ]

# Read in permutation results.
tajimasD_overall_norare_noncore$lower_p <- NA
tajimasD_overall_norare_noncore$upper_p <- NA

for (row_i in 1:nrow(tajimasD_overall_norare_noncore)) {
  
  sp <- tajimasD_overall_norare_noncore[row_i, "Species"]
  gene <- tajimasD_overall_norare_noncore[row_i, "Gene"]
  mean_tajimasD_val <- tajimasD_overall_norare_noncore[row_i, "mean_D"]
  
  if (is.na(mean_tajimasD_val)) { next }
  
  permutation_result_file <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/gene_haplotypes_tajimasD/permutation_outputs/",
                                   sp, "/", gene, "_rep_mean_results.tsv", sep = "")
  
  permutation_results <- read.table(permutation_result_file, header = TRUE, sep = "\t", row.names = 1)
  
  permutation_results <- permutation_results[which(! is.na(permutation_results$mean_D)), ]
  
  if (nrow(permutation_results) < 100) { next }
  
  tajimasD_overall_norare_noncore[row_i, "lower_p"] <- (length(which(permutation_results$mean_D <= mean_tajimasD_val)) + 1) / (nrow(permutation_results) + 1)
  
  tajimasD_overall_norare_noncore[row_i, "upper_p"] <- (length(which(permutation_results$mean_D >= mean_tajimasD_val)) + 1) / (nrow(permutation_results) + 1)
  
}

# Clear bias towards lower Tajima's D compared with permuted expectation

ggplot(data = tajimasD_overall_norare_noncore, aes(x = lower_p)) +
  geom_histogram() +
  theme_bw() +
  ylab("Number of genes") +
  xlab("p-value") +
  ggtitle("Prop. replicates where Tajima's D <= Observed")
  

all_corr_p <- p.adjust(c(tajimasD_overall_norare_noncore$lower_p, tajimasD_overall_norare_noncore$upper_p), "fdr")

tajimasD_overall_norare_noncore$lower_p_corr <- all_corr_p[1:nrow(tajimasD_overall_norare_noncore)]
tajimasD_overall_norare_noncore$upper_p_corr <- all_corr_p[(nrow(tajimasD_overall_norare_noncore) + 1):(nrow(tajimasD_overall_norare_noncore) * 2)]

# Sanity check:
which(tajimasD_overall_norare_noncore$lower_p_corr < 0.15 & tajimasD_overall_norare_noncore$upper_p_corr < 0.15)

all_lower_tajimasD_sig_genes <- tajimasD_overall_norare_noncore[which(tajimasD_overall_norare_noncore$lower_p_corr < 0.15), "Gene"]
all_upper_tajimasD_sig_genes <- tajimasD_overall_norare_noncore[which(tajimasD_overall_norare_noncore$upper_p_corr < 0.15), "Gene"]

all_panaroo <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/combined_panaroo_annot.rds")

all_lower_tajimasD_sig_genes_annot <- all_panaroo[all_lower_tajimasD_sig_genes, "Annotation", drop = FALSE]

all_upper_tajimasD_sig_genes_annot <- all_panaroo[all_upper_tajimasD_sig_genes, "Annotation", drop = FALSE]

sort(table(all_lower_tajimasD_sig_genes_annot$Annotation), decreasing = FALSE)
sort(table(all_upper_tajimasD_sig_genes_annot$Annotation), decreasing = FALSE)

# Species breakdown
table(tajimasD_overall_norare_noncore[which(tajimasD_overall_norare_noncore$lower_p_corr < 0.15), "Species"])
table(tajimasD_overall_norare_noncore[which(tajimasD_overall_norare_noncore$upper_p_corr < 0.15), "Species"])

# Show example

# Chose Bartonella_apis_aatA (Aspartate/prephenate aminotransferase)

lower_example_permutation_out <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/gene_haplotypes_tajimasD/permutation_outputs/Bartonella_apis/Bartonella_apis_aatA_rep_mean_results.tsv",
                                            header = TRUE, sep = "\t", row.names = 1)

ggplot(data = lower_example_permutation_out, aes(mean_D)) +
  geom_histogram(bins=100) +
  theme_bw() +
  ylab("Number of replicates") +
  xlab("Mean Tajima's D across samples") +
  ggtitle("Bartonella apis; aatA;\nAspartate/prephenate aminotransferase") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = tajimasD_overall_norare_noncore[which(tajimasD_overall_norare_noncore$Gene == "Bartonella_apis_aatA"), "mean_D"],
             linetype = 2,
             color = "red",
             size = 1)


# Run functional enrichment for all COG categories.
category_to_gene <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/COG_category_to_ortholog.rds")

background_sp_breakdown <- table(as.character(tajimasD_overall_norare_noncore[, "Species"]))

sp_upper_enrichments_raw <- list()

upper_sig_sp_breakdown <- table(as.character(tajimasD_overall_norare_noncore[which(tajimasD_overall_norare_noncore$upper_p_corr < 0.15), "Species"]))
upper_sig_sp_prevalent <- names(upper_sig_sp_breakdown)[which(upper_sig_sp_breakdown >= 10)]

for (sp in upper_sig_sp_prevalent) {
  
  sp_upper_sig_genes <- grep(sp, all_upper_tajimasD_sig_genes, value = TRUE)
  
  sp_other_tested_genes <- grep(sp, tajimasD_overall_norare_noncore$Gene, value = TRUE)
  sp_other_tested_genes <- sp_other_tested_genes[which(! sp_other_tested_genes %in% sp_upper_sig_genes)]
  
  sp_enrichment_tab <- identify_enriched_categories(genes = sp_upper_sig_genes, background = sp_other_tested_genes, category_to_gene_map = category_to_gene)
  
  # Only consider categories that have gene hits in at 10 genes (in either focal or background set).
  
  total_gene_hits <- sp_enrichment_tab$genes_num_category + sp_enrichment_tab$background_num_category
  
  sp_enrichment_tab <- sp_enrichment_tab[which(total_gene_hits >= 10), ]
  
  sp_enrichment_tab$fdr <- p.adjust(sp_enrichment_tab$p, "BH")
  
  sp_enrichment_tab$Species <- sp
  
  sp_upper_enrichments_raw[[sp]] <- sp_enrichment_tab

}

sp_upper_enrichments <- do.call(rbind, sp_upper_enrichments_raw)

sp_upper_enrichments$Description <- COG_category_descrip[sp_upper_enrichments$category, 1]

sp_upper_enrichments[which(sp_upper_enrichments$fdr < 0.3), ]




sp_lower_enrichments_raw <- list()

lower_sig_sp_breakdown <- table(as.character(tajimasD_overall_norare_noncore[which(tajimasD_overall_norare_noncore$lower_p_corr < 0.15), "Species"]))
lower_sig_sp_prevalent <- names(lower_sig_sp_breakdown)[which(lower_sig_sp_breakdown >= 10)]

for (sp in lower_sig_sp_prevalent) {
  
  sp_lower_sig_genes <- grep(sp, all_lower_tajimasD_sig_genes, value = TRUE)
  
  sp_other_tested_genes <- grep(sp, tajimasD_overall_norare_noncore$Gene, value = TRUE)
  sp_other_tested_genes <- sp_other_tested_genes[which(! sp_other_tested_genes %in% sp_lower_sig_genes)]
  
  sp_enrichment_tab <- identify_enriched_categories(genes = sp_lower_sig_genes, background = sp_other_tested_genes, category_to_gene_map = category_to_gene)
  
  # Only consider categories that have gene hits in at 10 genes (in either focal or background set).
  
  total_gene_hits <- sp_enrichment_tab$genes_num_category + sp_enrichment_tab$background_num_category
  
  sp_enrichment_tab <- sp_enrichment_tab[which(total_gene_hits >= 10), ]
  
  sp_enrichment_tab$fdr <- p.adjust(sp_enrichment_tab$p, "BH")
  
  sp_enrichment_tab$Species <- sp
  
  sp_lower_enrichments_raw[[sp]] <- sp_enrichment_tab
}

sp_lower_enrichments <- do.call(rbind, sp_lower_enrichments_raw)

sp_lower_enrichments$Description <- COG_category_descrip[sp_lower_enrichments$category, 1]

sp_lower_enrichments[which(sp_lower_enrichments$fdr < 0.15), c("Species", "category", "Description", "OR", "fdr")]


# Species     category                                                  Description        OR          fdr
# Bartonella_apis.25                               Bartonella_apis any_category                                                         <NA> 1.7306159 3.966390e-02
# Bifidobacterium_asteroides.16         Bifidobacterium_asteroides            G                        Carbohydrate transport and metabolism 1.5767129 5.922265e-02
# Bifidobacterium_asteroides.25         Bifidobacterium_asteroides any_category                                                         <NA> 1.9146629 3.623588e-02
# Gilliamella_apicola.21                       Gilliamella_apicola            P                       Inorganic ion transport and metabolism 7.0276243 3.059840e-02
# Gilliamella_apicola.25                       Gilliamella_apicola any_category                                                         <NA> 2.4529915 1.104346e-03
# Lactobacillus_apis.3                          Lactobacillus_apis            K                                                Transcription 0.5602979 9.192426e-02
# Lactobacillus_apis.13                         Lactobacillus_apis            O Posttranslational modification, protein turnover, chaperones 0.4339431 1.320794e-01
# Lactobacillus_apis.15                         Lactobacillus_apis            C                             Energy production and conversion 2.1995825 1.320794e-01
# Lactobacillus_apis.19                         Lactobacillus_apis            H                            Coenzyme transport and metabolism 4.6828753 2.375214e-03
# Lactobacillus_apis.25                         Lactobacillus_apis any_category                                                         <NA> 1.9946475 1.161386e-02
# Lactobacillus_helsingborgensis.25 Lactobacillus_helsingborgensis any_category                                                         <NA> 2.4153923 7.666534e-05
# Lactobacillus_melliventris.3          Lactobacillus_melliventris            K                                                Transcription 0.5228365 1.801319e-02
# Lactobacillus_melliventris.25         Lactobacillus_melliventris any_category                                                         <NA> 1.7002406 1.801319e-02




# % of genes identified as lower and upper dN/dS outliers (of the set tested with # alleles above the cut-off)
# Exclude species with < 100 tested genes for the purpose of this analysis.
sp2exclude <- names(background_sp_breakdown)[which(background_sp_breakdown < 100)]

upper_sig_sp_breakdown_percent <- (upper_sig_sp_breakdown / background_sp_breakdown[names(upper_sig_sp_breakdown)]) * 100
upper_sig_sp_breakdown_percent <- upper_sig_sp_breakdown_percent[which(! names(upper_sig_sp_breakdown_percent) %in% sp2exclude)]


lower_sig_sp_breakdown_percent <- (lower_sig_sp_breakdown / background_sp_breakdown[names(lower_sig_sp_breakdown)]) * 100
lower_sig_sp_breakdown_percent <- lower_sig_sp_breakdown_percent[which(! names(lower_sig_sp_breakdown_percent) %in% sp2exclude)]


sig_sp_breakdown_percent <- data.frame(species = names(upper_sig_sp_breakdown_percent),
                                       upper = as.numeric(upper_sig_sp_breakdown_percent),
                                       lower = as.numeric(lower_sig_sp_breakdown_percent[names(upper_sig_sp_breakdown_percent)]))

sig_sp_breakdown_percent$species <- gsub("_", " ", sig_sp_breakdown_percent$species)

ggplot(sig_sp_breakdown_percent, aes(x=upper, y=lower, label = species)) +
  geom_point(colour = "grey", alpha = 0.7, size = 3) +
  geom_text(colour="black", size = 3) +
  theme_bw() +
  xlab("% applicable genes with significantly high Tajima's D") +
  ylab("% applicable genes with significantly low Tajima's D") +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 60))

cor.test(sig_sp_breakdown_percent$upper, sig_sp_breakdown_percent$lower, method = "spearman")



# Look at how putatively highly mobile genes compare to others based on this test.

putatively_highly_mobile_genes <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/putative_mobile_genes_by_RVadj.tsv",
                                             header = FALSE, stringsAsFactors = FALSE)$V1

highly_mobile_tajimasD_output <- tajimasD_overall_norare_noncore[which(tajimasD_overall_norare_noncore$Gene %in% putatively_highly_mobile_genes), ]

hist(highly_mobile_tajimasD_output$lower_p_corr)
hist(highly_mobile_tajimasD_output$upper_p_corr)

highly_mobile_tajimasD_output[which(highly_mobile_tajimasD_output$upper_p_corr < 0.15), "Gene"]

highly_mobile_tajimasD_output[which(highly_mobile_tajimasD_output$lower_p_corr < 0.15 | highly_mobile_tajimasD_output$upper_p_corr < 0.15), "Gene"]

# 40/98

combined_panaroo[highly_mobile_tajimasD_output[which(highly_mobile_tajimasD_output$upper_p_corr < 0.15), "Gene"],]
combined_panaroo[highly_mobile_tajimasD_output[which(highly_mobile_tajimasD_output$lower_p_corr < 0.15), "Gene"], "Annotation"]