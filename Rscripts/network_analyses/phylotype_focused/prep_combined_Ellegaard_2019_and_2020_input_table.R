### Remove samples with insufficient coverage per genus from the input file for pandora.

rm(list = ls(all.names = TRUE))

pandora_output <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina/pandora_multisample.matrix",
                             header = TRUE, sep = "\t", row.names = 1)
rownames(pandora_output) <- gsub("\\.fa$", "", rownames(pandora_output))

# Identify core genes

panaroo_summary_and_core_genes <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/combined_panaroo_and_core_genes.rds")

# Table representing % of core genes called as present per phylotype per sample.
core_gene_coverage <- data.frame(matrix(NA, ncol = 5, nrow = 74))
colnames(core_gene_coverage) <- c("Bifidobacterium", "Firm4", "Firm5", "Gilliamella", "Snodgrassella")
rownames(core_gene_coverage) <- colnames(pandora_output)

for (phylotype in colnames(core_gene_coverage)) {
  
  phylotype_core_genes <- grep(phylotype, panaroo_summary_and_core_genes$core_genes, value = TRUE)
  
  for (samplename in rownames(core_gene_coverage)) {
    core_gene_coverage[samplename, phylotype] <- length(which(pandora_output[phylotype_core_genes, samplename] == 1)) / length(phylotype_core_genes) * 100
  }
}

max_coverage_per_sample <- apply(core_gene_coverage, 1, max, na.rm = TRUE)

# Remove three obvious outliers:
plot(max_coverage_per_sample, colSums(pandora_output))

colnames(pandora_output)[which(colSums(pandora_output) < 2000)]

pandora_output_filt <- pandora_output[, -which(colSums(pandora_output) < 2000)]

write.table(x = pandora_output_filt,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina/pandora_multisample.filt.matrix",
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)


pandora_output_filt_noncore <- pandora_output_filt[-which(rownames(pandora_output_filt) %in% panaroo_summary_and_core_genes$core_genes), ]

write.table(x = pandora_output_filt_noncore,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina/pandora_multisample.filt.noncore.matrix",
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

