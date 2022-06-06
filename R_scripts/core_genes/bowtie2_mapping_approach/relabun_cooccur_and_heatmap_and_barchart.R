### Test for correlations in species relative abundance across samples. 

rm(list = ls(all.names = TRUE))

library(Hmisc)
library(ComplexHeatmap)
library(ggplot2)
library(reshape2)

species_mean_depth_rel <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/summary/species_rel_abun.tsv.gz",
                                     header = TRUE, row.names = 1, sep = "\t")

species_mean_depth_rel_cor <- rcorr(as.matrix(species_mean_depth_rel), type = "spearman")

species_mean_depth_rel_cor$P[upper.tri(species_mean_depth_rel_cor$P)] <- NA
species_mean_depth_rel_cor$FDR <- p.adjust(species_mean_depth_rel_cor$P, "BH")

species_mean_depth_rel_cor$r_sig <- species_mean_depth_rel_cor$r
species_mean_depth_rel_cor$r_sig[! species_mean_depth_rel_cor$FDR < 0.25] <- NA
species_mean_depth_rel_cor$r_sig[is.na(species_mean_depth_rel_cor$FDR)] <- NA


Heatmap(matrix = species_mean_depth_rel_cor$r_sig, cluster_rows = FALSE, cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(! is.na(species_mean_depth_rel_cor$r_sig[i, j] > 0))
          grid.text(sprintf("%.1f", species_mean_depth_rel_cor$r_sig[i, j]), x, y, gp = gpar(fontsize = 10))
        })

par(mfrow = c(1, 2))
plot(species_mean_depth_rel$Bartonella_apis, species_mean_depth_rel$Bifidobacterium_asteroides)
plot(species_mean_depth_rel$Bartonella_apis, species_mean_depth_rel$Gilliamella_apicola)

plot(species_mean_depth_rel$Snodgrassella_alvi, species_mean_depth_rel$Lactobacillus_kimbladii)
plot(species_mean_depth_rel$Snodgrassella_alvi, species_mean_depth_rel$Lactobacillus_kullabergensis)

# Stacked barchart.

quick_hclust <- hclust(dist(species_mean_depth_rel))
quick_hclust$labels[quick_hclust$order]

species_relabun_long <- melt(as.matrix(species_mean_depth_rel))
species_relabun_long <- species_relabun_long[-which(species_relabun_long$value == 0), ]

colnames(species_relabun_long) <- c("Sample", "Species", "Relabun")

species_relabun_long$Sample <- factor(species_relabun_long$Sample, levels = quick_hclust$labels[quick_hclust$order])

plot_colours <- c("#d9407b",
                 "#5bc05d",
                 "#9b58cb",
                 "#b7b436",
                 "#5a6bcc",
                 "#6f9935",
                 "#cd56b4",
                 "#51945f",
                 "#d4463d",
                 "#4bc0b4",
                 "#d98c36",
                 "#5e93cd",
                 "#ba6342",
                 "#ab7bc0",
                 "#796c29",
                 "#e3849d",
                 "#c5aa68",
                 "#9f4868")

ggplot(data = species_relabun_long, aes(y = Sample, x = Relabun, fill = Species)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = plot_colours) +
  xlab("Relative abundance (%)") +
  ylab("Sample") +
  theme_bw()

