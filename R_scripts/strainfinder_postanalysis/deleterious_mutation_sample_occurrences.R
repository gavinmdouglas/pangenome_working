rm(list = ls(all.names = TRUE))

allele_sample_occurrences <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/allele_sample_occurrences.tsv",
                                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(allele_sample_occurrences) <- allele_sample_occurrences$Allele

allele_indels <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/gene_haplotype_deleterious_mutations/indels_out.tsv",
                            header = FALSE, sep = "\t")
colnames(allele_indels) <- c("gene", "haplotype", "nonframe_insert", "nonframe_del", "frame_insert", "frame_del")
rownames(allele_indels) <- paste(allele_indels$gene, as.character(allele_indels$haplotype), sep = "|||")

del_point_mutations <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/gene_haplotype_deleterious_mutations/start_stop_out.tsv",
                                  header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(del_point_mutations) <- c('name', 'allele_raw', 'canonical_start_codon_missing', 'canonical_stop_codon_missing', 'start_position',
                                   'leading_percent_truncated', 'premature_stop_position', 'expected_stop_position',
                                   'trailing_percent_truncated')


del_point_mutations$allele_raw <- gsub("^haplotype", "", del_point_mutations$allele_raw)
del_point_mutations$allele_raw <- gsub("_.*$", "", del_point_mutations$allele_raw)
rownames(del_point_mutations) <- paste(del_point_mutations$name, del_point_mutations$allele_raw, sep = "|||")



# Categorize alleles by deleterious mutation content.
indel_and_del_categories <- list()
indel_and_del_categories[["nonframe_insert"]] <- rownames(allele_indels)[which(allele_indels$nonframe_insert > 0)]
indel_and_del_categories[["nonframe_del_alleles"]] <- rownames(allele_indels)[which(allele_indels$nonframe_del > 0)]
indel_and_del_categories[["frame_insert_alleles"]] <- rownames(allele_indels)[which(allele_indels$frame_insert > 0)]
indel_and_del_categories[["frame_del_alleles"]] <- rownames(allele_indels)[which(allele_indels$frame_del > 0)]

indel_and_del_categories[["lost_start_alleles"]] <- rownames(del_point_mutations)[which(del_point_mutations$canonical_start_codon_missing == "True")]
indel_and_del_categories[["lost_stop_alleles"]] <- rownames(del_point_mutations)[which(del_point_mutations$canonical_stop_codon_missing == "True")]
indel_and_del_categories[["premature_stop_alleles"]] <- rownames(del_point_mutations)[which(! is.na(del_point_mutations$premature_stop_position))]


category_subsets <- list()

for (category in names(indel_and_del_categories)) {
  category_subsets[[category]] <- allele_sample_occurrences[indel_and_del_categories[[category]], ]
  category_subsets[[category]]$allele_category <- category
}

allele_sample_occurrences_w_categories <- do.call(rbind, category_subsets)

allele_sample_occurrences_no_category <- allele_sample_occurrences[which(! rownames(allele_sample_occurrences) %in% allele_sample_occurrences_w_categories), ]
allele_sample_occurrences_no_category$allele_category <- "Other"

allele_sample_occurrences_w_categories <- rbind(allele_sample_occurrences_w_categories, allele_sample_occurrences_no_category)

allele_sample_occurrences_w_categories$allele_category[which(allele_sample_occurrences_w_categories$allele_category %in% c("frame_insert_alleles", "frame_del_alleles"))] <- "Frameshift"
allele_sample_occurrences_w_categories$allele_category[which(allele_sample_occurrences_w_categories$allele_category %in% c("nonframe_insert", "nonframe_del_alleles"))] <- "Non-frameshift"

allele_sample_occurrences_w_categories$allele_category[which(allele_sample_occurrences_w_categories$allele_category == "Other")] <- "All other alleles"

# For now just consider INDELs (and background)
allele_sample_occurrences_w_categories <- allele_sample_occurrences_w_categories[which(allele_sample_occurrences_w_categories$allele_category %in% c("Frameshift", "Non-frameshift", "All other alleles")), ]

allele_sample_occurrences_w_categories$allele_category <- factor(allele_sample_occurrences_w_categories$allele_category,
                                                                 levels = c("Frameshift", "Non-frameshift", "All other alleles"))

ggplot(data = allele_sample_occurrences_w_categories, aes(y = allele_category, x = Sample_occurrences)) +
  geom_quasirandom(groupOnX = FALSE) +
  geom_boxplot(fill = "grey", alpha = 0.6) +
  theme_bw() +
  xlab("Sample occurrences") +
  facet_wrap(Species ~ ., scales = "free_x") +
  ylab("")
