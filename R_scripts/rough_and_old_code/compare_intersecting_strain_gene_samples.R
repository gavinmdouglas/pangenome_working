rm(list = ls(all.names = TRUE))

total_num_samples <- 71

core_haplotype_abun <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_genes_output/strainfinder_pandora_struct_abun.rds")

all_gene_samples <- list.files(path = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/prepped_input/",
                               pattern = "samples",
                               full.names = TRUE)

all_gene_samples_no_path <- list.files(path = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/prepped_input/",
                                       pattern = "samples",
                                       full.names = FALSE)
all_gene_samples <- all_gene_samples[grep("no.struct", all_gene_samples)]

all_gene_samples_no_path <- all_gene_samples_no_path[grep("no.struct", all_gene_samples_no_path)]

all_genes <- gsub("_samples.txt", "", all_gene_samples_no_path)

sample_concordance <- data.frame(matrix(NA, nrow = length(all_genes), ncol = 7))
colnames(sample_concordance) <- c("phylotype", "gene", "num_gene_samples", "num_strain_samples", "minor_count", "exp_intersect", "observed_intersect")

rownames(sample_concordance) <- all_genes
sample_concordance$gene <- all_genes


for (i in 1:length(all_genes)) {

  gene_samples <- read.table(file = all_gene_samples[i], stringsAsFactors = FALSE)$V1
  
  gene <- all_genes[i]
  
  phylotype <- strsplit(gene, split = "_")[[1]][1]
  phylotype_samples <- rownames(core_haplotype_abun[[phylotype]])
  
  sample_concordance[gene, "phylotype"] <- phylotype
  
  sample_concordance[gene, c("num_gene_samples", "num_strain_samples")] <- c(length(gene_samples), length(phylotype_samples))
  
  sample_concordance[gene, "minor_count"] <- min(c(length(gene_samples), length(phylotype_samples)))
  
  sample_concordance[gene, "exp_intersect"] <- (length(gene_samples) / total_num_samples) * (length(phylotype_samples) / total_num_samples) * total_num_samples
  
  sample_concordance[gene, "observed_intersect"] <- length(which(gene_samples %in% phylotype_samples))
}

library(ggplot2)
library(reshape2)

sample_concordance <- sample_concordance[-which(sample_concordance$phylotype == "Firm4"), ]

sample_concordance$observed_intersect_by_minor_count <- sample_concordance$observed_intersect / sample_concordance$minor_count


ggplot(data = sample_concordance, aes(y = observed_intersect, x = exp_intersect, colour = phylotype)) +
  geom_point(size=3) +
  geom_abline(intercept = 0, slope = 1, linetype=2)

hist_out <- ggplot(data = sample_concordance, aes(x = observed_intersect_by_minor_count, fill = phylotype)) +
  geom_histogram() +
  facet_grid(. ~ phylotype)

scatterplot_out <- ggplot(data = sample_concordance, aes(x = observed_intersect_by_minor_count, y = exp_intersect, colour = phylotype)) +
  geom_point(size = 3) +
  facet_grid(. ~ phylotype)

plot_grid(hist_out, scatterplot_out, nrow = 2)
 