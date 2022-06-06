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

sp_dnds_overall_raw <- list()

dnds_working_dir <- "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/gene_haplotypes_dNdS/"

for (sp in species_present) {
 
  dnds_summary_files <- list.files(path = paste(dnds_working_dir, sp, sep = "/"), full.names = TRUE, pattern = "_dnds.tsv$") 
  dnds_summary_files <- dnds_summary_files[grep("pairwise_haplotype_dnds.tsv$", dnds_summary_files, invert = TRUE)]
  
  sp_dnds_overall <- data.frame(matrix(NA, nrow = length(dnds_summary_files), ncol = 5))
  colnames(sp_dnds_overall) <- c("Species", "Gene", "all_alleles", "within_dnds", "between_dnds")
  sp_genes <- basename(gsub("_dnds.tsv", "", dnds_summary_files))
  rownames(sp_dnds_overall) <- sp_genes
  
  sp_dnds_overall$Species <- sp
  sp_dnds_overall$Gene <- sp_genes
  
  for (i in 1:length(dnds_summary_files)) {
    raw_mean_dnds <- read.table(dnds_summary_files[i], header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
    
    dnds_samples_present <- rownames(raw_mean_dnds)[which(rownames(raw_mean_dnds) %in% all_samples)]
    
    sp_dnds_overall[i, "all_alleles"] <- raw_mean_dnds["all_haplotypes", "within_dnds"]
    sp_dnds_overall[i, "within_dnds"] <- raw_mean_dnds["all_samples", "within_dnds"]
    sp_dnds_overall[i, "between_dnds"] <- raw_mean_dnds["all_samples", "between_dnds"]
    
  }
  
  sp_dnds_overall_raw[[sp]] <- sp_dnds_overall

}

dnds_overall <- do.call(rbind, sp_dnds_overall_raw)

dnds_overall$Species <- factor(dnds_overall$Species)
dnds_overall$num_alleles <- num_alleles[dnds_overall$Gene, ]

dnds_overall_long <- melt(dnds_overall, ids = c("Species", "Gene"), variable.name = "Type", value.name = "dN_dS")

dnds_overall$ratio <- dnds_overall$within_dnds / dnds_overall$between_dnds

# Exclude genes with fewer than 6 alleles (as a permutation test based on that many allelic combos is almost certainly not useful).
# Chosen because there are 720 different permutations of 6 elements, whereas there are only 120 permutations of 5 elements.
dnds_overall_norare <- dnds_overall[which(dnds_overall$num_alleles >= 6), ]

# Read in permutation results.
dnds_overall_norare$lower_p <- NA
dnds_overall_norare$upper_p <- NA

for (row_i in 1:nrow(dnds_overall_norare)) {
  
  sp <- dnds_overall_norare[row_i, "Species"]
  gene <- dnds_overall_norare[row_i, "Gene"]
  within_dnds_val <- dnds_overall_norare[row_i, "within_dnds"]
  
  if (is.na(within_dnds_val)) { next }
  
  permutation_result_file <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/gene_haplotypes_dNdS/permutation_outputs/",
                                   sp, "/", gene, "_rep_mean_dnds.tsv", sep = "")
  
  permutation_results <- read.table(permutation_result_file, header = TRUE, sep = "\t", row.names = 1)
  
  permutation_results <- permutation_results[which(! is.na(permutation_results$within_dnds)), ]
  
  if (nrow(permutation_results) < 100) { next }
  
  dnds_overall_norare[row_i, "lower_p"] <- (length(which(permutation_results$within_dnds <= within_dnds_val)) + 1) / (nrow(permutation_results) + 1)
  
  dnds_overall_norare[row_i, "upper_p"] <- (length(which(permutation_results$within_dnds >= within_dnds_val)) + 1) / (nrow(permutation_results) + 1)
  
}

all_corr_p <- p.adjust(c(dnds_overall_norare$lower_p, dnds_overall_norare$upper_p), "fdr")

dnds_overall_norare$lower_p_corr <- all_corr_p[1:nrow(dnds_overall_norare)]
dnds_overall_norare$upper_p_corr <- all_corr_p[(nrow(dnds_overall_norare) + 1):(nrow(dnds_overall_norare) * 2)]

# Sanity check:
which(dnds_overall_norare$lower_p_corr < 0.15 & dnds_overall_norare$upper_p_corr < 0.15)

all_lower_dnds_sig_genes <- dnds_overall_norare[which(dnds_overall_norare$lower_p_corr < 0.15), "Gene"]
all_upper_dnds_sig_genes <- dnds_overall_norare[which(dnds_overall_norare$upper_p_corr < 0.15), "Gene"]

all_panaroo <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/combined_panaroo_annot.rds")

all_lower_dnds_sig_genes_annot <- all_panaroo[all_lower_dnds_sig_genes, "Annotation", drop = FALSE]

all_upper_dnds_sig_genes_annot <- all_panaroo[all_upper_dnds_sig_genes, "Annotation", drop = FALSE]

sort(table(all_lower_dnds_sig_genes_annot$Annotation), decreasing = FALSE)
sort(table(all_upper_dnds_sig_genes_annot$Annotation), decreasing = FALSE)

# Species breakdown
table(dnds_overall_norare[which(dnds_overall_norare$lower_p_corr < 0.15), "Species"])
table(dnds_overall_norare[which(dnds_overall_norare$upper_p_corr < 0.15), "Species"])

# Show examples

all_upper_dnds_sig_genes_annot[which(all_upper_dnds_sig_genes_annot$Annotation == "FAD:protein FMN transferase"), , drop = FALSE]
# Chose Snodgrassella_alvi_apbE 

upper_example_permutation_out <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/gene_haplotypes_dNdS/permutation_outputs/Snodgrassella_alvi/Snodgrassella_alvi_apbE_rep_mean_dnds.tsv",
                                            header = TRUE, sep = "\t", row.names = 1)

dnds_overall[which(dnds_overall$Gene == "Snodgrassella_alvi_apbE"), ]

ggplot(data = upper_example_permutation_out, aes(within_dnds)) +
  geom_histogram(bins=100) +
  theme_bw() +
  ylab("Number of replicates") +
  xlab("Mean within-sample dN/dS") +
  ggtitle("Snodgrassella alvi; apbE;\nFAD:protein FMN transferase") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0.105, 0.175)) +
  geom_vline(xintercept = dnds_overall[which(dnds_overall$Gene == "Snodgrassella_alvi_apbE"), "within_dnds"],
             linetype = 2,
             color = "red",
             size = 1)

all_lower_dnds_sig_genes_annot[which(all_lower_dnds_sig_genes_annot$Annotation == "PTS system beta-glucoside-specific EIIBCA component"), , drop = FALSE]
# Chose Lactobacillus_apis_bglF_6 

lower_example_permutation_out <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/gene_haplotypes_dNdS/permutation_outputs/Lactobacillus_apis/Lactobacillus_apis_bglF_6_rep_mean_dnds.tsv",
                                            header = TRUE, sep = "\t", row.names = 1)

dnds_overall[which(dnds_overall$Gene == "Lactobacillus_apis_bglF_6"), ]

ggplot(data = lower_example_permutation_out, aes(within_dnds)) +
  geom_histogram(bins=100) +
  theme_bw() +
  ylab("Number of replicates") +
  xlab("Mean within-sample dN/dS") +
  ggtitle("Lactobacillus apis; bglF_6;\nPTS system beta-glucoside-specific EIIBCA component") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = dnds_overall[which(dnds_overall$Gene == "Lactobacillus_apis_bglF_6"), "within_dnds"],
             linetype = 2,
             color = "red",
             size = 1)


# Run functional enrichment for all COG categories.
category_to_gene <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/COG_category_to_ortholog.rds")

background_sp_breakdown <- table(as.character(dnds_overall_norare[, "Species"]))

sp_upper_enrichments_raw <- list()

upper_sig_sp_breakdown <- table(as.character(dnds_overall_norare[which(dnds_overall_norare$upper_p_corr < 0.15), "Species"]))
upper_sig_sp_prevalent <- names(upper_sig_sp_breakdown)[which(upper_sig_sp_breakdown >= 10)]

for (sp in upper_sig_sp_prevalent) {
  
  sp_upper_sig_genes <- grep(sp, all_upper_dnds_sig_genes, value = TRUE)
  
  sp_other_tested_genes <- grep(sp, dnds_overall_norare$Gene, value = TRUE)
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

sp_upper_enrichments[which(sp_upper_enrichments$fdr < 0.15), ]

# category genes_num_category genes_num_other background_num_category background_num_other         OR           p       fdr
# Bifidobacterium_asteroides.6         V                  4              25                      71                 1575  3.5492958 0.037637631 0.2634634
# Bifidobacterium_asteroides.12        U                  2              27                      10                 1636 12.1185185 0.017161194 0.2106851
# Bifidobacterium_asteroides.22        Q                  2              27                      11                 1635 11.0101010 0.020065247 0.2106851
# Bombilactobacillus_mellis.20         I                  8              72                      38                 1160  3.3918129 0.006153493 0.1292234
# Gilliamella_apicola.5                D                  6             120                      31                 1985  3.2016129 0.018953879 0.2427507
# Gilliamella_apicola.23               R                  4             122                     180                 1836  0.3344262 0.021108755 0.2427507
# Species                                                   Description
# Bifidobacterium_asteroides.6  Bifidobacterium_asteroides                                            Defense mechanisms
# Bifidobacterium_asteroides.12 Bifidobacterium_asteroides Intracellular trafficking, secretion, and vesicular transport
# Bifidobacterium_asteroides.22 Bifidobacterium_asteroides  Secondary metabolites biosynthesis, transport and catabolism
# Bombilactobacillus_mellis.20   Bombilactobacillus_mellis                                Lipid transport and metabolism
# Gilliamella_apicola.5                Gilliamella_apicola    Cell cycle control, cell division, chromosome partitioning
# Gilliamella_apicola.23               Gilliamella_apicola                              General function prediction only




sp_lower_enrichments_raw <- list()

lower_sig_sp_breakdown <- table(as.character(dnds_overall_norare[which(dnds_overall_norare$lower_p_corr < 0.15), "Species"]))
lower_sig_sp_prevalent <- names(lower_sig_sp_breakdown)[which(lower_sig_sp_breakdown >= 10)]

for (sp in lower_sig_sp_prevalent) {
  
  sp_lower_sig_genes <- grep(sp, all_lower_dnds_sig_genes, value = TRUE)
  
  sp_other_tested_genes <- grep(sp, dnds_overall_norare$Gene, value = TRUE)
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

sp_lower_enrichments[which(sp_lower_enrichments$fdr < 0.15), ]


# category genes_num_category genes_num_other background_num_category background_num_other       OR           p        fdr
# Bartonella_apis.8                        M                  9              40                     134                 1628 2.733582 0.012356294 0.25948218
# Gilliamella_apicola.20                   I                  7              63                      55                 2017 4.074747 0.003448092 0.07930612
# Lactobacillus_helsingborgensis.15        C                  6              43                      45                 1291 4.003101 0.007674822 0.14582161
# Species                            Description
# Bartonella_apis.8                                Bartonella_apis Cell wall/membrane/envelope biogenesis
# Gilliamella_apicola.20                       Gilliamella_apicola         Lipid transport and metabolism
# Lactobacillus_helsingborgensis.15 Lactobacillus_helsingborgensis       Energy production and conversion



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

sig_sp_breakdown_percent$species[which(sig_sp_breakdown_percent$species == "Lactobacillus_helsingborgensis")] <- "Lactobacillus\nhelsingborgensis"
sig_sp_breakdown_percent$species[which(sig_sp_breakdown_percent$species == "Gilliamella_apis")] <- "Gilliamella\napis"
sig_sp_breakdown_percent$species[which(sig_sp_breakdown_percent$species == "Bifidobacterium_coryneforme_indicum")] <- "Bifidobacterium_coryneforme/indicum"
sig_sp_breakdown_percent$species <- gsub("_", " ", sig_sp_breakdown_percent$species)

ggplot(sig_sp_breakdown_percent, aes(x=upper, y=lower, label = species)) +
  geom_point(colour = "grey", alpha = 0.7, size = 3) +
  geom_text(colour="black", size = 3) +
  theme_bw() +
  xlab("% applicable genes with significantly high dN/dS") +
  ylab("% applicable genes with significantly low dN/dS")

cor.test(sig_sp_breakdown_percent$upper, sig_sp_breakdown_percent$lower, method = "spearman")


sig_sp_breakdown_percent$num_strains <- sapply(names(upper_sig_sp_breakdown_percent), function(x) { ncol(strain_abun[[x]]) })
sig_sp_breakdown_percent$num_samples <- sapply(names(upper_sig_sp_breakdown_percent), function(x) { nrow(strain_abun[[x]]) })


sig_sp_breakdown_percent$summed_upper_lower <- sig_sp_breakdown_percent$upper + sig_sp_breakdown_percent$lower

ggplot(sig_sp_breakdown_percent, aes(x=num_strains, y=summed_upper_lower, label = species)) +
  geom_point(colour = "grey", alpha = 0.7, size = 3) +
  geom_text(colour="black", size = 3) +
  theme_bw() +
  xlab("# strains") +
  ylab("% applicable genes with significantly different dN/dS")

ggplot(sig_sp_breakdown_percent, aes(x=num_samples, y=summed_upper_lower, label = species)) +
  geom_point(colour = "grey", alpha = 0.7, size = 3) +
  theme_bw() +
  xlab("# strains") +
  ylab("% applicable genes with significantly different dN/dS")






# Look at how putatively highly mobile genes compare to others based on this test.

putatively_highly_mobile_genes <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/putative_mobile_genes_by_RVadj.tsv",
                                             header = FALSE, stringsAsFactors = FALSE)$V1

highly_mobile_dnds_output <- dnds_overall_norare[which(dnds_overall_norare$Gene %in% putatively_highly_mobile_genes), ]

hist(highly_mobile_dnds_output$lower_p_corr)
hist(highly_mobile_dnds_output$upper_p_corr)

highly_mobile_dnds_output[which(highly_mobile_dnds_output$lower_p_corr < 0.15 | highly_mobile_dnds_output$upper_p_corr < 0.15), "Gene"]

# Only 6/98

