rm(list = ls(all.names = TRUE))

# Dig into RVadj outputs.

library(ggplot2)
library(ggbeeswarm)

strain_abun <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/strain_abun.rds")

all_genes_abun <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/all_genes_OTU.rds")


RVadj_out <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/strain_vs_gene_abun_RVadj.rds")

length(which(is.na(RVadj_out$RVadj)))

RVadj_out <- RVadj_out[which(! is.na(RVadj_out$RVadj)), ]

# First, summarize across all species and genes (irrespective of how many samples they're in / if a permutation test was done).
ggplot(data = RVadj_out, aes(y = species, x = RVadj)) +
  geom_quasirandom(groupOnX = FALSE, col = "grey") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  theme_bw()


# Restricted only to genes with at least 20 intersecting samples (which is the subset I restricted the permutation tests to)

RVadj_out_20intersect <- RVadj_out[which(RVadj_out$num_intersecting >= 20), ]

RVadj_out_20intersect$permuted_RVadj_p_set <- "Non sig."
RVadj_out_20intersect$permuted_RVadj_p_set[which(RVadj_out_20intersect$permuted_RVadj_p < 0.05)] <- "Sig."

# Added in # strains and # gene haplotypes per comparison, to see how that relates to the results.

RVadj_out_20intersect$strain_num <- sapply(RVadj_out_20intersect$species, function(x) { ncol(strain_abun[[x]]) })
RVadj_out_20intersect$allele_num <- sapply(RVadj_out_20intersect$gene,
                                           function(x) {
                                             x_split <- strsplit(x = x, "_")[[1]]
                                             sp <- paste(x_split[1], x_split[2], sep = "_")
                                             
                                             if (sp == "Bifidobacterium_coryneforme") { sp <- "Bifidobacterium_coryneforme_indicum" }
                                             
                                             ncol(all_genes_abun[[sp]][[x]])
                                           })


ggplot(data = RVadj_out_20intersect, aes(y = RVadj, x = allele_num, col = permuted_RVadj_p_set)) +
  geom_point() +
  facet_wrap(. ~ species) +
  theme_bw() +
  xlab("# alleles") +
  ylab("RVadj") +
  ggtitle("All genes - restricted to cases with > 20 intersecting samples")
  


RVadj_out_20intersect_10alleles <- RVadj_out_20intersect[which(RVadj_out_20intersect$allele_num >= 10), ]


RVadj_out_20intersect_10alleles$permuted_RVadj_p_col <- "black"
RVadj_out_20intersect_10alleles$permuted_RVadj_p_col[which(RVadj_out_20intersect_10alleles$permuted_RVadj_p_set == "Sig.")] <- "red"

ggplot(data = RVadj_out_20intersect_10alleles, aes(y = species, x = RVadj)) +
  geom_quasirandom(groupOnX = FALSE, col = RVadj_out_20intersect_10alleles$permuted_RVadj_p_col) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  theme_bw() +
  ggtitle("All genes - restricted to cases with > 20 intersecting samples\nand >= 10 alleles")




# Breakdown of proportion sig/non-sig. genes per species

uniq_species <- sort(unique(RVadj_out_20intersect_10alleles$species))
prop_sig <- data.frame(matrix(NA, ncol = 2, nrow = length(uniq_species)))
rownames(prop_sig) <- uniq_species
colnames(prop_sig) <- c("Sig.", "Non-sig.")

for (sp in uniq_species) {
  sp_subset <- RVadj_out_20intersect_10alleles[which(RVadj_out_20intersect_10alleles$species == sp), ]
  prop_sig[sp, "Sig."] <- length(which(sp_subset$permuted_RVadj_p < 0.05))
  prop_sig[sp, "Non-sig."] <- length(which(sp_subset$permuted_RVadj_p >= 0.05))
}

prop_sig$prop_sig <- prop_sig$Sig. / (prop_sig$`Non-sig.` + prop_sig$Sig.)


prop_sig$num_samples_w_species <- sapply(rownames(prop_sig), function(x) { nrow(strain_abun[[x]]) })

plot(prop_sig$num_samples_w_species, prop_sig$prop_sig, xlim = c(0, 74), ylim = c(0, 1),
     xlab = "Species prevalence (based on strains across samples)", ylab = "Prop. of genes that were significant")

rownames(prop_sig)[which(prop_sig$prop_sig > 0.9)]




# Restrict to species that reached that prevalence plateau.

prevalent_species <- rownames(prop_sig)[which(prop_sig$prop_sig > 0.9)]

RVadj_out_20intersect_10alleles_prevalent <- RVadj_out_20intersect_10alleles[which(RVadj_out_20intersect_10alleles$species %in% prevalent_species), ]
table(RVadj_out_20intersect_10alleles_prevalent[which(RVadj_out_20intersect_10alleles_prevalent$permuted_RVadj_p_set == "Non sig."), "species"])


# Save this list of putatively mobile genes to a file.

putative_mobile_genes <- RVadj_out_20intersect_10alleles_prevalent[which(RVadj_out_20intersect_10alleles_prevalent$permuted_RVadj_p_set == "Non sig."), "gene"]

write.table(x = putative_mobile_genes, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/putative_mobile_genes_by_RVadj.tsv",
            col.names = FALSE, row.names = FALSE, quote = FALSE)


# Run functional enrichment, by COG categories.

source("/home/gdouglas/scripts/pangenome_working/R_scripts/functions.R")

category_to_gene <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/COG_category_to_ortholog.rds")


all_sig_genes <- RVadj_out_20intersect_10alleles_prevalent[which(RVadj_out_20intersect_10alleles_prevalent$permuted_RVadj_p_set == "Sig."), "gene"]

sp_enrichments_raw <- list()

for (sp in prevalent_species) {

  sp_putative_mobile_genes <- grep(sp, putative_mobile_genes, value = TRUE)

  sp_sig_genes <- grep(sp, all_sig_genes, value = TRUE)

  sp_enrichment_tab <- identify_enriched_categories(genes = sp_putative_mobile_genes, background = sp_sig_genes, category_to_gene_map = category_to_gene)

  # Only consider categories that have gene hits in at 10 genes (in either focal or background set).
  
  total_gene_hits <- sp_enrichment_tab$genes_num_category + sp_enrichment_tab$background_num_category
  
  sp_enrichment_tab <- sp_enrichment_tab[which(total_gene_hits >= 10), ]

  sp_enrichment_tab$fdr <- p.adjust(sp_enrichment_tab$p, "BH")

  sp_enrichment_tab$Species <- sp
    
  sp_enrichments_raw[[sp]] <- sp_enrichment_tab
}

sp_enrichments <- do.call(rbind, sp_enrichments_raw)


sig_enrichments <- sp_enrichments[which(sp_enrichments$fdr < 0.05), ]
rownames(sig_enrichments) <- NULL



# Get table detailing functional annotations for each of these putatively mobile genes.

panaroo_out <- list()

for (sp in prevalent_species) {
   
  panaroo_out[[sp]] <- read.table(paste("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/", sp, "/gene_presence_absence.csv.gz", sep = ""),
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
  
  panaroo_out[[sp]] <- panaroo_out[[sp]][, c(1:3)]
  
  rownames(panaroo_out[[sp]]) <- paste(sp, panaroo_out[[sp]]$Gene, sep = "_")
  
}

panaroo_out <- do.call(rbind, panaroo_out)
rownames(panaroo_out) <- gsub("^.*\\.", "", rownames(panaroo_out))

panaroo_out_putative_mobile <- panaroo_out[putative_mobile_genes, ]