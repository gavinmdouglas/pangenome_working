rm(list = ls(all.names = TRUE))

# Explore if the distribution of different alleles per gene departs from a
# basic model where every strain increases the probability that a different allele will be present proportionally.

library(ggplot2)
library(ggbeeswarm)
library(cowplot)

source("/home/gdouglas/scripts/pangenome_working/R_scripts/functions.R")



strain_abun <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/strain_abun.rds")

num_strains_per_species <- sapply(strain_abun, nrow)

species_w_10_strains <- names(num_strains_per_species)[which(num_strains_per_species >= 10)]

all_genes_abun <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/all_genes_OTU.rds")

# Exclude core genes
panaroo_only_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/panaroo_only_potential_core.rds")
checkm_panaroo_passed_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/checkm_panaroo_passed_potential_core.rds")
core_genes <- c(panaroo_only_core, checkm_panaroo_passed_core)
all_core_genes <- do.call(c, core_genes)
names(all_core_genes) <- NULL

for (sp in names(all_genes_abun)) {
 noncore_genes <- names(all_genes_abun[[sp]])[which(! names(all_genes_abun[[sp]]) %in% all_core_genes)]
 all_genes_abun[[sp]] <- all_genes_abun[[sp]][noncore_genes]
}

all_sp_allele_info <- list()

for (sp in species_w_10_strains) {
  
  test_strain <- strain_abun[[sp]]
  num_strains <- ncol(test_strain)
  num_strains_per_sample <- rowSums(test_strain > 0)
  
  
  all_num_alleles <- as.numeric()
  sample_prevalence <- as.numeric()
  mean_num_alleles_when_present <- as.numeric()
  max_num_alleles_when_present <- as.numeric()
  mean_ratio_num_alleles_strains <- as.numeric()

  full_p_out <- as.numeric()
  cor_rho_out <- as.numeric()
  cor_p_out <- as.numeric()
  
  cor_present_rho_out <- as.numeric()
  cor_present_p_out <- as.numeric()
  
  for (g in names(all_genes_abun[[sp]])) {
  
    test_gene <- all_genes_abun[[sp]][[g]]
  
    num_alleles <- ncol(test_gene)
    
    all_num_alleles <- c(all_num_alleles, num_alleles)
    
    num_alleles_per_sample <- rowSums(test_gene > 0)
  
    mean_num_alleles_when_present <- c(mean_num_alleles_when_present, mean(num_alleles_per_sample))
    max_num_alleles_when_present <- c(max_num_alleles_when_present, max(num_alleles_per_sample))
    
    test_gene <- test_gene[which(rownames(test_gene) %in% rownames(test_strain)), , drop = FALSE]
    test_gene <- test_gene[, which(colSums(test_gene) > 0), drop = FALSE]
    num_alleles_per_sample <- rowSums(test_gene > 0)
    
    if (nrow(test_gene) == 0 | num_alleles <= 2) {
      sample_prevalence <- c(sample_prevalence, 0)
      full_p_out <- c(full_p_out, NA)
      
      cor_rho_out <- c(cor_rho_out, NA)
      cor_p_out <- c(cor_p_out, NA)
      
      cor_present_rho_out <- c(cor_present_rho_out, NA)
      cor_present_p_out <- c(cor_present_p_out, NA)
      
      mean_ratio_num_alleles_strains <- c(mean_ratio_num_alleles_strains, NA)
      
      next
    }
    
    strain_samples_missing_alleles <- names(num_strains_per_sample)[which(! names(num_strains_per_sample) %in% names(num_alleles_per_sample))]
    
    if (length(strain_samples_missing_alleles) > 0) {
      missing_vec <- rep(0, length(strain_samples_missing_alleles))
      names(missing_vec) <- strain_samples_missing_alleles
      num_alleles_per_strain_samples <- c(num_alleles_per_sample, missing_vec)
      num_alleles_per_strain_samples <- num_alleles_per_strain_samples[names(num_strains_per_sample)]
    }
    
    sample_prevalence <- c(sample_prevalence, length(which(num_alleles_per_strain_samples > 0)))
    
    mean_ratio_num_alleles_strains <- c(mean_ratio_num_alleles_strains, mean(num_alleles_per_sample / num_strains_per_sample[names(num_alleles_per_sample)]))
    
    # Multinomial test based on all samples (including those where the gene is absent).
    exp_prob <- sum(num_alleles_per_strain_samples) / sum(num_strains_per_sample)
    exp_lik <- exp_prob * num_strains_per_sample
    exp_prob_dist <- exp_lik / sum(exp_lik)
    full_p <- XNomial::xmonte(obs = num_alleles_per_strain_samples, expr = exp_prob_dist, ntrials = 1000)$pLLR
    full_p_out <- c(full_p_out, full_p)
  
    cor_out <- cor.test(num_alleles_per_strain_samples, num_strains_per_sample, method = "spearman")
    cor_rho_out <- c(cor_rho_out, cor_out$estimate)
    cor_p_out <- c(cor_p_out, cor_out$p.value)
    
    if (length(num_alleles_per_sample) >= 10) { 
      cor_present_out <- cor.test(num_alleles_per_sample, num_strains_per_sample[names(num_alleles_per_sample)], method = "spearman")
      cor_present_rho_out <- c(cor_present_rho_out, cor_present_out$estimate)
      cor_present_p_out <- c(cor_present_p_out, cor_present_out$p.value)
    } else{
      cor_present_rho_out <- c(cor_present_rho_out, NA)
      cor_present_p_out <- c(cor_present_p_out, NA)
    }
  }
  
  sp_allele_info <- data.frame(Species = sp,
                               Gene = names(all_genes_abun[[sp]]),
                               num_alleles = all_num_alleles,
                               sample_prev = sample_prevalence,
                               mean_num_alleles_when_present = mean_num_alleles_when_present,
                               max_num_alleles_when_present = max_num_alleles_when_present,
                               mean_ratio_num_alleles_strains_when_present = mean_ratio_num_alleles_strains,
                               full_p = full_p_out,
                               rho = cor_rho_out,
                               cor_p =  cor_p_out,
                               rho_present = cor_present_rho_out,
                               cor_present_p = cor_present_p_out)
  
  all_sp_allele_info[[sp]] <- sp_allele_info                          

}

allele_info <- do.call(rbind, all_sp_allele_info)

# Set min p-value to be 0.0009
allele_info[which(allele_info$full_p == 0), "full_p"] <- 0.0009
allele_info$full_fdr <- p.adjust(allele_info$full_p, "BH")


# For the correlation, ignore all genes that are present in < 10 samples with strains, or that have < 10 unique alleles, as well.
allele_info$cor_p[which(allele_info$sample_prev < 10 | allele_info$num_alleles < 10)] <- NA
allele_info$rho[which(allele_info$sample_prev < 10 | allele_info$num_alleles < 10)] <- NA
allele_info$cor_fdr <- p.adjust(allele_info$cor_p, "BH")


allele_info$cor_present_p[which(allele_info$sample_prev < 10 | allele_info$num_alleles < 10)] <- NA
allele_info$rho_present[which(allele_info$sample_prev < 10 | allele_info$num_alleles < 10)] <- NA
allele_info$cor_present_fdr <- p.adjust(allele_info$cor_present_p, "BH")


# Breakdown of spearman correlations.
spearman_results_raw <- allele_info[, c("Species", "Gene", "rho", "cor_p", "cor_fdr")]
spearman_results_raw$type <- "All samples"
colnames(spearman_results_raw) <- c("Species", "Gene", "Rho", "p", "FDR", "Type")

# Get mean per species
species_mean_rho <- aggregate(Rho ~ Species, data = spearman_results_raw, FUN = mean)
species_mean_rho <- species_mean_rho[order(species_mean_rho$Rho), ]

spearman_results_present <- allele_info[, c("Species", "Gene", "rho_present", "cor_present_p", "cor_present_fdr")]
spearman_results_present$type <- "Samples with gene present only"
colnames(spearman_results_present) <- c("Species", "Gene", "Rho", "p", "FDR", "Type")

spearman_results <- rbind(spearman_results_raw, spearman_results_present)

spearman_results <- spearman_results[-which(is.na(spearman_results$Rho) | is.na(spearman_results$p)), ]

spearman_results$sig_col <- "red"
spearman_results$sig_col[which(spearman_results$p >= 0.05)] <- "grey"

spearman_results$Species <- factor(spearman_results$Species, levels = species_mean_rho$Species)

ggplot(data = spearman_results, aes(x = Rho, y = Species, colour = p)) +
  geom_quasirandom(colour = spearman_results$sig_col, groupOnX = FALSE, dodge.width=.75) +
  geom_boxplot(fill = "grey", outlier.shape = NA, alpha = 0.7) +
  xlab("Spearman's rho (# alleles vs # strains)") +
  theme_bw() +
  facet_grid(. ~ Type)



# Species breakdown
sp_gene_counts <- table(allele_info$Species)
sp_sig_gene_counts <- table(allele_info[which(allele_info$full_fdr < 0.05), "Species"])
species_sig_breakdown <- data.frame(species = names(sp_gene_counts),
                                    gene_count = as.integer(sp_gene_counts),
                                    sig_gene_count = as.integer(sp_sig_gene_counts[names(sp_gene_counts)]))

species_sig_breakdown$percent_sig <- (species_sig_breakdown$sig_gene_count / species_sig_breakdown$gene_count) * 100

species_sig_breakdown <- species_sig_breakdown[order(species_sig_breakdown$percent_sig), ]

species_sig_breakdown$species <- factor(species_sig_breakdown$species, levels = rev(species_sig_breakdown$species))

ggplot(data = species_sig_breakdown, aes(x = percent_sig, y = species)) +
  geom_bar(stat="identity") +
  theme_bw() +
  geom_text(label = as.character(round(species_sig_breakdown$percent_sig, 1)), nudge_x = 3) +
  ylab("") +
  xlab("% tested non-core genes sig.") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 105))


# COG category enrichment for these significant genes.
COG_category_descrip <- read.table("/data1/gdouglas/db/COG_definitions/COG_category_descrip.tsv",
                                   header = FALSE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
colnames(COG_category_descrip) <- "COG_category"

COG_category_descrip["any_category", "COG_category"] <- "any_category"


category_to_gene <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/COG_category_to_ortholog.rds")

sp_enrichments_raw <- list()

for (sp in species_sig_breakdown$species) {
  
  sp_sig_genes <- allele_info[which(allele_info$Species == sp & allele_info$full_fdr < 0.05), "Gene"]
  
  sp_nonsig_genes <- allele_info[which(allele_info$Species == sp & allele_info$full_fdr >= 0.05), "Gene"]

  
  sp_enrichment_tab <- identify_enriched_categories(genes = sp_sig_genes, background = sp_nonsig_genes, category_to_gene_map = category_to_gene)
  
  # Only consider categories that have gene hits in at 10 genes (in either focal or background set).
  
  total_gene_hits <- sp_enrichment_tab$genes_num_category + sp_enrichment_tab$background_num_category
  
  sp_enrichment_tab <- sp_enrichment_tab[which(total_gene_hits >= 10), ]
  
  sp_enrichment_tab$fdr <- p.adjust(sp_enrichment_tab$p, "BH")
  
  sp_enrichment_tab$Species <- sp
  
  sp_enrichments_raw[[sp]] <- sp_enrichment_tab
}

sp_enrichments <- do.call(rbind, sp_enrichments_raw)

sp_enrichments$Description <- COG_category_descrip[sp_enrichments$category, 1]

sp_enrichments[which(sp_enrichments$fdr < 0.15), c("OR", "fdr", "Description")]

# OR          fdr                                                   Description
# Bartonella_apis.6                 3.00711913 7.734762e-02                                            Defense mechanisms
# Bartonella_apis.25                0.65048810 7.734762e-02                                                  any_category
# Lactobacillus_helsingborgensis.6  2.52404372 1.201609e-01                                            Defense mechanisms
# Lactobacillus_helsingborgensis.7  0.55108179 1.463059e-01                                Signal transduction mechanisms
# Lactobacillus_helsingborgensis.17 0.61600631 1.463059e-01                           Amino acid transport and metabolism
# Lactobacillus_helsingborgensis.21 0.56732833 1.467639e-01                        Inorganic ion transport and metabolism
# Lactobacillus_helsingborgensis.25 0.53458935 5.360365e-02                                                  any_category
# Bombilactobacillus_mellis.25      1.54166667 1.037990e-01                                                  any_category
# Lactobacillus_apis.1              0.43902439 8.938812e-02               Translation, ribosomal structure and biogenesis
# Lactobacillus_apis.3              1.82510524 1.212539e-01                                                 Transcription
# Lactobacillus_apis.5              0.32891216 8.938812e-02    Cell cycle control, cell division, chromosome partitioning
# Lactobacillus_apis.6              2.49110629 1.381133e-01                                            Defense mechanisms
# Lactobacillus_apis.8              0.52097130 4.693381e-02                        Cell wall/membrane/envelope biogenesis
# Lactobacillus_apis.17             0.48937008 8.363555e-03                           Amino acid transport and metabolism
# Lactobacillus_apis.25             0.34527972 7.979402e-05                                                  any_category
# Gilliamella_apicola.1             0.08746095 7.969825e-02               Translation, ribosomal structure and biogenesis
# Gilliamella_apicola.12            0.16508152 7.969825e-02 Intracellular trafficking, secretion, and vesicular transport
# Bifidobacterium_asteroides.1      0.13369963 4.640489e-02               Translation, ribosomal structure and biogenesis
# Bifidobacterium_asteroides.6      0.32476471 4.640489e-02                                            Defense mechanisms
# Bifidobacterium_asteroides.15     0.08780488 1.289449e-02                              Energy production and conversion
# Bifidobacterium_asteroides.21     0.35131744 1.072693e-01                        Inorganic ion transport and metabolism
# Lactobacillus_melliventris.3      2.65918896 2.477672e-02                                                 Transcription
# Lactobacillus_melliventris.7      0.52377115 1.101993e-01                                Signal transduction mechanisms
# Lactobacillus_melliventris.15     0.42504550 6.320720e-02                              Energy production and conversion
# Lactobacillus_melliventris.17     0.55978992 7.613834e-02                           Amino acid transport and metabolism
# Lactobacillus_melliventris.21     0.47560976 6.180113e-02                        Inorganic ion transport and metabolism
# Lactobacillus_melliventris.25     0.34743802 2.890669e-04                                                  any_category



# Show examples
allele_info[which(allele_info$Gene == "Bartonella_apis_arsC"), ]
# Arsenate reductase

sp <- "Bartonella_apis"
g <- "Bartonella_apis_arsC"
example_strains_counts <- rowSums(strain_abun[[sp]] > 0)
example_allele_counts <- rowSums(all_genes_abun[[sp]][[g]] > 0)
intersecting_samples <- names(example_allele_counts)[which(names(example_allele_counts) %in% names(example_strains_counts))]
example_allele_counts <- example_allele_counts[intersecting_samples]

strain_samples_missing_alleles <- names(example_strains_counts)[which(! names(example_strains_counts) %in% names(example_allele_counts))]

if (length(strain_samples_missing_alleles) > 0) {
  missing_vec <- rep(0, length(strain_samples_missing_alleles))
  names(missing_vec) <- strain_samples_missing_alleles
  example_allele_counts <- c(example_allele_counts, missing_vec)
  example_allele_counts <- example_allele_counts[names(example_strains_counts)]
}

example_counts = data.frame(strains=example_strains_counts, alleles=example_allele_counts)

ex1 <- ggplot(data = example_counts, aes(x = strains, y = alleles)) +
  geom_point(position=position_jitter(h=0.15,w=0.15), pch = 1, fill = NA, size = 5) +
  theme_bw() +
  xlab("# strains") +
  ylab("# alleles") +
  ggtitle("Bartonella apis - arsC alleles") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(0, 20) +
  ylim(0, 10)




allele_info[which(allele_info$Gene == "Bifidobacterium_asteroides_arnC"), ]
# Catalyzes the transfer of 4-deoxy-4-formamido-L-arabinose from UDP to undecaprenyl phosphate.
# The modified arabinose is attached to lipid A and is required for resistance to polymyxin and cationic antimicrobial peptides.

sp <- "Bifidobacterium_asteroides"
g <- "Bifidobacterium_asteroides_arnC"
example_strains_counts <- rowSums(strain_abun[[sp]] > 0)
example_allele_counts <- rowSums(all_genes_abun[[sp]][[g]] > 0)
intersecting_samples <- names(example_allele_counts)[which(names(example_allele_counts) %in% names(example_strains_counts))]
example_allele_counts <- example_allele_counts[intersecting_samples]

strain_samples_missing_alleles <- names(example_strains_counts)[which(! names(example_strains_counts) %in% names(example_allele_counts))]

if (length(strain_samples_missing_alleles) > 0) {
  missing_vec <- rep(0, length(strain_samples_missing_alleles))
  names(missing_vec) <- strain_samples_missing_alleles
  example_allele_counts <- c(example_allele_counts, missing_vec)
  example_allele_counts <- example_allele_counts[names(example_strains_counts)]
}

example_counts = data.frame(strains=example_strains_counts, alleles=example_allele_counts)

ex2 <- ggplot(data = example_counts, aes(x = strains, y = alleles)) +
  geom_point(position=position_jitter(h=0.15,w=0.15), pch = 1, fill = NA, size = 5) +
  theme_bw() +
  xlab("# strains") +
  ylab("# alleles") +
  ggtitle("Bifidobacterium asteroides - arnC alleles") +
  theme(plot.title = element_text(hjust = 0.5))

plot_grid(ex1, ex2)




# Breakdown of mean ration of num alleles vs num strains (sanity check)

ggplot(data = allele_info, aes(y = Species, x = mean_ratio_num_alleles_strains_when_present)) +
  geom_quasirandom(groupOnX = FALSE, col = "grey") +
  geom_boxplot(outlier.shape = NA, col = "black", fill = "grey", alpha = 0.3) +
  theme_bw() +
  ylab("") +
  xlab("Mean (# alleles / # strains)") +
  scale_y_discrete(limits=rev) +
  geom_vline(xintercept = 1, lwd = 1, lty = 2, col = "red")

##### Unused code, used for data exploration #####

# # Quick look into whether there are examples of genes that are widely prevalent,
# # and for which there are multiple alleles, but the number of alleles per sample is low
# # Limit to cases with at least 5 alleles and in 10 samples
# 
# allele_info_atleast5_and_10samples <- allele_info[which(allele_info$num_alleles >= 5 & allele_info$sample_prev >= 10), ]
# 
# allele_info_atleast5_and_10samples$max_prop <- allele_info_atleast5_and_10samples$max_num_alleles_when_present / allele_info_atleast5_and_10samples$num_alleles
# allele_info_atleast5_and_10samples$mean_prop <- allele_info_atleast5_and_10samples$mean_num_alleles_when_present / allele_info_atleast5_and_10samples$num_alleles
# 
# plot(allele_info_atleast5_and_10samples$num_alleles, allele_info_atleast5_and_10samples$mean_prop)
# plot(allele_info_atleast5_and_10samples$mean_ratio_num_alleles_strains_when_present, allele_info_atleast5_and_10samples$mean_prop)
# plot(allele_info_atleast5_and_10samples$mean_num_alleles_when_present, allele_info_atleast5_and_10samples$mean_ratio_num_alleles_strains_when_present)
# 
# allele_info_atleast5[which(allele_info_atleast5_and_10samples$mean_prop < 0.3 & allele_info_atleast5_and_10samples$mean_ratio_num_alleles_strains_when_present < 0.5), "Gene"]
#                            
#                            
# # annot (exploring)
# 
# combined_panaroo <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/combined_panaroo_annot.rds")
# 
# table(combined_panaroo[allele_info_atleast5[which(allele_info_atleast5_and_10samples$mean_prop < 0.3 & allele_info_atleast5_and_10samples$mean_ratio_num_alleles_strains_when_present < 0.5), "Gene"], "Annotation"])
# 
# test_genes <- grep("Gilliamella", allele_info_atleast5[which(allele_info_atleast5_and_10samples$mean_prop < 0.3 & allele_info_atleast5_and_10samples$mean_ratio_num_alleles_strains_when_present < 0.5), "Gene"], value = TRUE)
# table(combined_panaroo[test_genes, "Annotation"])
# 
# combined_panaroo["Lactobacillus_helsingborgensis_benM_3~~~catM", "Annotation"]
