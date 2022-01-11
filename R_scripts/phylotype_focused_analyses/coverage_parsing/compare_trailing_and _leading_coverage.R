rm(list = ls(all.names = TRUE))

# Compare how depth changes at the beginning and end of genes as you move towards centre.
# This is based on tables that contain the average over all genes in a sample.
# The metric is the depth per site divided by the median depth of all sites with non-zero depth.

library(cowplot)
library(ggplot2)

setwd("honey_bee_pangenome/data/Ellegaard_panaroo_pangenome_mapping/mapped_genes_edge_median_ratios/")

# Process the Gilliamella data
Gilliamella_ratio_files <- list.files("stringent_comp_map/Gilliamella_per_bp_median_ratio", full.names = TRUE, pattern = ".gz")

Gilliamella_starting_file <- Gilliamella_ratio_files[1]
Gilliamella_ratio_files <- Gilliamella_ratio_files[-1]

Gilliamella_ratios <- read.table(Gilliamella_starting_file, header = TRUE, sep = "\t")
Gilliamella_starting_sample <- gsub("Gilliamella_per_bp_median_ratio/", "", Gilliamella_starting_file)
Gilliamella_starting_sample <- gsub(".Gilliamella.merged.nonparalog.per_site.edge_ratio.tsv.gz", "", Gilliamella_starting_sample)
Gilliamella_ratios$sample <- Gilliamella_starting_sample

for (ratio_file in Gilliamella_ratio_files) {
 
  Gilliamella_ratio_tab <- read.table(ratio_file, header = TRUE, sep = "\t")
  tmp_Gilliamella_sample <- gsub("Gilliamella_per_bp_median_ratio/", "", ratio_file)
  tmp_Gilliamella_sample <- gsub(".Gilliamella.merged.nonparalog.per_site.edge_ratio.tsv.gz", "", tmp_Gilliamella_sample)
  Gilliamella_ratio_tab$sample <- tmp_Gilliamella_sample

  Gilliamella_ratios <- rbind(Gilliamella_ratios, Gilliamella_ratio_tab)
  
}

Gilliamella_ratios_leading <- Gilliamella_ratios[which(Gilliamella_ratios$gene_position > 0), ]
Gilliamella_ratios_trailing <- Gilliamella_ratios[which(Gilliamella_ratios$gene_position < 0), ]

Gilliamella_ratios_leading_plot <- ggplot(data = Gilliamella_ratios_leading, aes(x = gene_position, y = mean_ratio, group = gene_position)) +
                                          geom_boxplot(outlier.shape = NA) +
                                          ggtitle(expression(paste(italic("Gilliamella"), " gene-leading depth"))) +
                                          xlab("Position in gene (from start)") +
                                          ylab("Normalized depth") +
                                          theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
                                          coord_cartesian(ylim = c(0, 1.6))

Gilliamella_ratios_trailing_plot <- ggplot(data = Gilliamella_ratios_trailing, aes(x = gene_position, y = mean_ratio, group = gene_position)) +
                                           geom_boxplot(outlier.shape = NA) +
                                           ggtitle(expression(paste(italic("Gilliamella"), " gene-trailing depth"))) +
                                           xlab("Position in gene (from end)") +
                                           ylab("Normalized depth") +
                                           theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
                                           coord_cartesian(ylim = c(0, 1.6))


# Process the Snodgrassella data
Snodgrassella_ratio_files <- list.files("stringent_comp_map/Snodgrassella_per_bp_median_ratio", full.names = TRUE, pattern = ".gz")

Snodgrassella_starting_file <- Snodgrassella_ratio_files[1]
Snodgrassella_ratio_files <- Snodgrassella_ratio_files[-1]

Snodgrassella_ratios <- read.table(Snodgrassella_starting_file, header = TRUE, sep = "\t")
Snodgrassella_starting_sample <- gsub("Snodgrassella_per_bp_median_ratio/", "", Snodgrassella_starting_file)
Snodgrassella_starting_sample <- gsub(".Snodgrassella.merged.nonparalog.per_site.edge_ratio.tsv.gz", "", Snodgrassella_starting_sample)
Snodgrassella_ratios$sample <- Snodgrassella_starting_sample

for (ratio_file in Snodgrassella_ratio_files) {
  
  Snodgrassella_ratio_tab <- read.table(ratio_file, header = TRUE, sep = "\t")
  tmp_Snodgrassella_sample <- gsub("Snodgrassella_per_bp_median_ratio/", "", ratio_file)
  tmp_Snodgrassella_sample <- gsub(".Snodgrassella.merged.nonparalog.per_site.edge_ratio.tsv.gz", "", tmp_Snodgrassella_sample)
  Snodgrassella_ratio_tab$sample <- tmp_Snodgrassella_sample
  
  Snodgrassella_ratios <- rbind(Snodgrassella_ratios, Snodgrassella_ratio_tab)
  
}

Snodgrassella_ratios_leading <- Snodgrassella_ratios[which(Snodgrassella_ratios$gene_position > 0), ]
Snodgrassella_ratios_trailing <- Snodgrassella_ratios[which(Snodgrassella_ratios$gene_position < 0), ]

Snodgrassella_ratios_leading_plot <- ggplot(data = Snodgrassella_ratios_leading, aes(x = gene_position, y = mean_ratio, group = gene_position)) +
                                          geom_boxplot(outlier.shape = NA) +
                                          ggtitle(expression(paste(italic("Snodgrassella"), " gene-leading depth"))) +
                                          xlab("Position in gene (from start)") +
                                          ylab("Normalized depth") +
                                          theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
                                          coord_cartesian(ylim = c(0, 1.6))

Snodgrassella_ratios_trailing_plot <- ggplot(data = Snodgrassella_ratios_trailing, aes(x = gene_position, y = mean_ratio, group = gene_position)) +
                                           geom_boxplot(outlier.shape = NA) +
                                           ggtitle(expression(paste(italic("Snodgrassella"), " gene-trailing depth"))) +
                                           xlab("Position in gene (from end)") +
                                           ylab("Normalized depth") +
                                           theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
                                           coord_cartesian(ylim = c(0, 1.6))



# Combine plots
plot_grid(Gilliamella_ratios_leading_plot,
          Gilliamella_ratios_trailing_plot,
          Snodgrassella_ratios_leading_plot,
          Snodgrassella_ratios_trailing_plot,
          labels = c('a', 'b', 'c', 'd'),
          nrow = 2,
          ncol = 2)








# As above but for the simple comp-mapped data
rm(list = ls(all.names = TRUE))

# Process the Gilliamella data
Gilliamella_ratio_files <- list.files("best_hit_simple_comp_map//Gilliamella_per_bp_median_ratio", full.names = TRUE, pattern = ".gz")

Gilliamella_starting_file <- Gilliamella_ratio_files[1]
Gilliamella_ratio_files <- Gilliamella_ratio_files[-1]

Gilliamella_ratios <- read.table(Gilliamella_starting_file, header = TRUE, sep = "\t")
Gilliamella_starting_sample <- gsub("Gilliamella_per_bp_median_ratio/", "", Gilliamella_starting_file)
Gilliamella_starting_sample <- gsub(".Gilliamella.merged.nonparalog.per_site.edge_ratio.tsv.gz", "", Gilliamella_starting_sample)
Gilliamella_ratios$sample <- Gilliamella_starting_sample

for (ratio_file in Gilliamella_ratio_files) {
  
  Gilliamella_ratio_tab <- read.table(ratio_file, header = TRUE, sep = "\t")
  tmp_Gilliamella_sample <- gsub("Gilliamella_per_bp_median_ratio/", "", ratio_file)
  tmp_Gilliamella_sample <- gsub(".Gilliamella.merged.nonparalog.per_site.edge_ratio.tsv.gz", "", tmp_Gilliamella_sample)
  Gilliamella_ratio_tab$sample <- tmp_Gilliamella_sample
  
  Gilliamella_ratios <- rbind(Gilliamella_ratios, Gilliamella_ratio_tab)
  
}

Gilliamella_ratios_leading <- Gilliamella_ratios[which(Gilliamella_ratios$gene_position > 0), ]
Gilliamella_ratios_trailing <- Gilliamella_ratios[which(Gilliamella_ratios$gene_position < 0), ]

Gilliamella_ratios_leading_plot <- ggplot(data = Gilliamella_ratios_leading, aes(x = gene_position, y = mean_ratio, group = gene_position)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle(expression(paste(italic("Gilliamella"), " gene-leading depth"))) +
  xlab("Position in gene (from start)") +
  ylab("Normalized depth") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
  coord_cartesian(ylim = c(0, 1.6))

Gilliamella_ratios_trailing_plot <- ggplot(data = Gilliamella_ratios_trailing, aes(x = gene_position, y = mean_ratio, group = gene_position)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle(expression(paste(italic("Gilliamella"), " gene-trailing depth"))) +
  xlab("Position in gene (from end)") +
  ylab("Normalized depth") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
  coord_cartesian(ylim = c(0, 1.6))


# Process the Snodgrassella data
Snodgrassella_ratio_files <- list.files("best_hit_simple_comp_map/Snodgrassella_per_bp_median_ratio", full.names = TRUE, pattern = ".gz")

Snodgrassella_starting_file <- Snodgrassella_ratio_files[1]
Snodgrassella_ratio_files <- Snodgrassella_ratio_files[-1]

Snodgrassella_ratios <- read.table(Snodgrassella_starting_file, header = TRUE, sep = "\t")
Snodgrassella_starting_sample <- gsub("Snodgrassella_per_bp_median_ratio/", "", Snodgrassella_starting_file)
Snodgrassella_starting_sample <- gsub(".Snodgrassella.merged.nonparalog.per_site.edge_ratio.tsv.gz", "", Snodgrassella_starting_sample)
Snodgrassella_ratios$sample <- Snodgrassella_starting_sample

for (ratio_file in Snodgrassella_ratio_files) {
  
  Snodgrassella_ratio_tab <- read.table(ratio_file, header = TRUE, sep = "\t")
  tmp_Snodgrassella_sample <- gsub("Snodgrassella_per_bp_median_ratio/", "", ratio_file)
  tmp_Snodgrassella_sample <- gsub(".Snodgrassella.merged.nonparalog.per_site.edge_ratio.tsv.gz", "", tmp_Snodgrassella_sample)
  Snodgrassella_ratio_tab$sample <- tmp_Snodgrassella_sample
  
  Snodgrassella_ratios <- rbind(Snodgrassella_ratios, Snodgrassella_ratio_tab)
  
}

Snodgrassella_ratios_leading <- Snodgrassella_ratios[which(Snodgrassella_ratios$gene_position > 0), ]
Snodgrassella_ratios_trailing <- Snodgrassella_ratios[which(Snodgrassella_ratios$gene_position < 0), ]

Snodgrassella_ratios_leading_plot <- ggplot(data = Snodgrassella_ratios_leading, aes(x = gene_position, y = mean_ratio, group = gene_position)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle(expression(paste(italic("Snodgrassella"), " gene-leading depth"))) +
  xlab("Position in gene (from start)") +
  ylab("Normalized depth") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
  coord_cartesian(ylim = c(0, 1.6))

Snodgrassella_ratios_trailing_plot <- ggplot(data = Snodgrassella_ratios_trailing, aes(x = gene_position, y = mean_ratio, group = gene_position)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle(expression(paste(italic("Snodgrassella"), " gene-trailing depth"))) +
  xlab("Position in gene (from end)") +
  ylab("Normalized depth") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
  coord_cartesian(ylim = c(0, 1.6))



# Combine plots
plot_grid(Gilliamella_ratios_leading_plot,
          Gilliamella_ratios_trailing_plot,
          Snodgrassella_ratios_leading_plot,
          Snodgrassella_ratios_trailing_plot,
          labels = c('a', 'b', 'c', 'd'),
          nrow = 2,
          ncol = 2)

