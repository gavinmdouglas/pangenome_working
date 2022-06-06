rm(list = ls(all.names = TRUE))

library(ggplot2)

species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species_present.txt",
                      stringsAsFactors = FALSE)$V1

all_tested_genes <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/all_genes_prepped.txt",
                               stringsAsFactors = FALSE)$V1

allele_sample_occurrences <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/allele_sample_occurrences.tsv",
                                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(allele_sample_occurrences) <- allele_sample_occurrences$Allele

allele_indels <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/gene_haplotype_deleterious_mutations/indels_out.tsv",
                            header = FALSE, sep = "\t")
colnames(allele_indels) <- c("gene", "haplotype", "nonframe_insert", "nonframe_del", "frame_insert", "frame_del")

allele_indels_summed <- aggregate(. ~ gene, data = allele_indels[, -2], FUN =sum)

nonframe_insert_genes <- allele_indels_summed$gene[which(allele_indels_summed$nonframe_insert > 0)]
nonframe_del_genes <- allele_indels_summed$gene[which(allele_indels_summed$nonframe_del > 0)]
frame_insert_genes <- allele_indels_summed$gene[which(allele_indels_summed$frame_insert > 0)]
frame_del_genes <- allele_indels_summed$gene[which(allele_indels_summed$frame_del > 0)]

indel_prop_summary <- data.frame(matrix(NA, nrow = length(species), ncol = 6))
rownames(indel_prop_summary) <- species
colnames(indel_prop_summary) <- c("nonframe_insert", "nonframe_del", "nonframe_either",
                                  "frame_insert", "frame_del", "frame_either")

for (sp in species) {
  sp_genes <- grep(sp, all_tested_genes, value = TRUE)
  
  indel_prop_summary[sp, "nonframe_insert"] <- length(which(sp_genes %in% nonframe_insert_genes)) / length(sp_genes)
  indel_prop_summary[sp, "nonframe_del"] <- length(which(sp_genes %in% nonframe_del_genes)) / length(sp_genes)
  indel_prop_summary[sp, "nonframe_either"] <- length(which(sp_genes %in% nonframe_insert_genes | sp_genes %in% nonframe_del_genes)) / length(sp_genes)
  
  indel_prop_summary[sp, "frame_insert"] <- length(which(sp_genes %in% frame_insert_genes)) / length(sp_genes)
  indel_prop_summary[sp, "frame_del"] <- length(which(sp_genes %in% frame_del_genes)) / length(sp_genes)
  indel_prop_summary[sp, "frame_either"] <- length(which(sp_genes %in% frame_insert_genes | sp_genes %in% frame_del_genes)) / length(sp_genes)
  
}

indel_percent_summary <- indel_prop_summary * 100

indel_percent_summary$Species <- rownames(indel_percent_summary)

ggplot(data = indel_percent_summary, aes(x = nonframe_either, y = frame_either)) +
  geom_point(size = 8, alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, lwd = 1, lty = 2) +
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 20)) +
  theme_bw() +
  xlab("% genes with non-frameshift insertion/deletion variants") +
  ylab("% genes with frameshift insertion/deletion variants")



# Get SFS of genes with frameshifts vs those with non-frameshifts
allele_sample_occurrences$frameshift

rownames(allele_indels) <- paste(allele_indels$gene, as.character(allele_indels$haplotype), sep = "|||")

