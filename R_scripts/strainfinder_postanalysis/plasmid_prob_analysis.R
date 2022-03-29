# Investigate plasmid probability of genes.

# In particular, check if carbohydrate-related genes (COG category G) have higher plasmid probabilities compared with other genes.
# This is relevant as these genes were highly enriched in the "putatively highly mobile set"

rm(list = ls(all.names = TRUE))

source("/home/gdouglas/scripts/pangenome_working/R_scripts/functions.R")

library(ggplot2)
library(ggbeeswarm)
library(cowplot)

plasmid_prob <- read.table("/data1/gdouglas/projects/honey_bee/ref_genomes/plasclass_out/plasclass_output_ortholog_summary.tsv",
                           header = TRUE, sep = "\t", row.names = 1)
plasmid_prob$Species <- sapply(rownames(plasmid_prob), function(x) { paste(strsplit(x, "_")[[1]][c(1, 2)], collapse = "_") })
plasmid_prob$Species[which(plasmid_prob$Species == "Bifidobacterium_coryneforme")] <- "Bifidobacterium_coryneforme_indicum"
plasmid_prob$Species <- factor(plasmid_prob$Species)

mean_plasclass_by_species <- ggplot(data = plasmid_prob, aes(x = mean, y = Species)) +
                                    geom_quasirandom(groupOnX = FALSE, size = 0.5) +
                                    geom_boxplot(fill = "grey", outlier.shape = NA, alpha = 0.7, lwd = 1) +
                                    theme_bw() +
                                    scale_y_discrete(limits = rev(levels(plasmid_prob$Species))) +
                                    xlab("Mean plasclass probability per ortholog") +
                                    ggtitle("Mean plasmid probability") +
                                    ylab("") +
                                    theme(plot.title = element_text(hjust = 0.5))

max_plasclass_by_species <- ggplot(data = plasmid_prob, aes(x = max, y = Species)) +
                                    geom_quasirandom(groupOnX = FALSE, size = 0.5) +
                                    geom_boxplot(fill = "grey", outlier.shape = NA, alpha = 0.7, lwd = 1) +
                                    theme_bw() +
                                    scale_y_discrete(limits = rev(levels(plasmid_prob$Species))) +
                                    xlab("Max plasclass probability per ortholog") +
                                    ggtitle("Max plasmid probability") +
                                    ylab("") +
                                    theme(plot.title = element_text(hjust = 0.5))

plot_grid(mean_plasclass_by_species, max_plasclass_by_species)



# Compare plasclass probability by COG category.
gene_COG_categories <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/all_microbiota_COG.category_by_ortholog.tsv",
                                  header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
colnames(gene_COG_categories) <- "COG_category"

COG_category_descrip <- read.table("/data1/gdouglas/db/COG_definitions/COG_category_descrip.tsv",
                                   header = FALSE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
colnames(COG_category_descrip) <- "COG_category"

COG_category_descrip["-", "COG_category"] <- "No COG"

plasmid_prob_w_COG <- plasmid_prob
plasmid_prob_w_COG$COG_category <- "-"
genes_w_COGs <- rownames(plasmid_prob_w_COG)[which(rownames(plasmid_prob_w_COG) %in% rownames(gene_COG_categories))]
plasmid_prob_w_COG[genes_w_COGs, "COG_category"] <- gene_COG_categories[genes_w_COGs, "COG_category"]
plasmid_prob_w_COG <- expand_for_multi_COG_category(plasmid_prob_w_COG)

# Remove species / COG category combinations where there are fewer than 10 non-NA values.
for (species in unique(plasmid_prob_w_COG$Species)) {
  
  for (category in unique(plasmid_prob_w_COG$COG_category)) {
    
    plasmid_prob_w_COG_species <- plasmid_prob_w_COG[which(plasmid_prob_w_COG$Species == species & plasmid_prob_w_COG$COG_category == category), "mean"]
    
    if (length(plasmid_prob_w_COG_species) == 0) { next }
    
    if (length(which(! is.na(plasmid_prob_w_COG_species))) < 10) {
      plasmid_prob_w_COG <- plasmid_prob_w_COG[-which(plasmid_prob_w_COG$Species == species & plasmid_prob_w_COG$COG_category == category), ]
    }
    
  }
  
}


ggplot(data = plasmid_prob_w_COG, aes(x = COG_category, y = mean)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey") +
  theme_bw() +
  ylab("Mean plasclass probability per ortholog") +
  xlab("COG category") +
  facet_wrap(Species ~ .)


# Focused specifically on Bifidobacterium asteroides
plasmid_prob_w_COG_Bifidobacterium_asteroides <- plasmid_prob_w_COG

plasmid_prob_w_COG_Bifidobacterium_asteroides <- plasmid_prob_w_COG_Bifidobacterium_asteroides[which(plasmid_prob_w_COG_Bifidobacterium_asteroides$Species == "Bifidobacterium_asteroides"), ]

ggplot(data = plasmid_prob_w_COG_Bifidobacterium_asteroides, aes(x = COG_category, y = mean)) +
  geom_quasirandom(groupOnX = TRUE, size = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey") +
  theme_bw() +
  ylab("Mean plasclass probability per ortholog") +
  xlab("COG category") +
  facet_wrap(Species ~ .)



# Then specifically look at categories enriched for putatively highly mobile genes and compare genes in this set vs all others.

putatively_high_mobilility_genes <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/putative_mobile_genes_by_RVadj.tsv",
                                                header = FALSE, sep = "\t")$V1

Bifidobacterium_asteroides_G_genes <- plasmid_prob_w_COG_Bifidobacterium_asteroides[which(plasmid_prob_w_COG_Bifidobacterium_asteroides$COG_category == "G"), ]
Bifidobacterium_asteroides_G_genes$RVadj_agreement <- "Sig."

Bifidobacterium_asteroides_G_genes[which(rownames(Bifidobacterium_asteroides_G_genes) %in% putatively_high_mobilility_genes), "RVadj_agreement"] <- "Non-sig."

ggplot(data = Bifidobacterium_asteroides_G_genes, aes(x = RVadj_agreement, y = mean))+
  geom_quasirandom(groupOnX = TRUE) +
  theme_bw() +
  ylab("Mean plasmid probability") +
  xlab("RVadj agreement with strain abundances") +
  ggtitle("Bifidobacterium asteroides carbohydrate-related genes (G)") +
  theme(plot.title = element_text(hjust = 0.5))


wilcox.test(Bifidobacterium_asteroides_G_genes$mean[which(Bifidobacterium_asteroides_G_genes$RVadj_agreement == "Non-sig.")],
            Bifidobacterium_asteroides_G_genes$mean[which(Bifidobacterium_asteroides_G_genes$RVadj_agreement == "Sig.")])


Bifidobacterium_asteroides_noCOG_genes <- plasmid_prob_w_COG_Bifidobacterium_asteroides[which(plasmid_prob_w_COG_Bifidobacterium_asteroides$COG_category == "-"), ]
Bifidobacterium_asteroides_noCOG_genes$RVadj_agreement <- "Sig."

Bifidobacterium_asteroides_noCOG_genes[which(rownames(Bifidobacterium_asteroides_noCOG_genes) %in% putatively_high_mobilility_genes), "RVadj_agreement"] <- "Non-sig."

ggplot(data = Bifidobacterium_asteroides_noCOG_genes, aes(x = RVadj_agreement, y = mean)) +
  geom_quasirandom(groupOnX = TRUE) +
  theme_bw() +
  ylab("Mean plasmid probability") +
  xlab("RVadj agreement with strain abundances") +
  ggtitle("Bifidobacterium asteroides non-COG-associated genes") +
  theme(plot.title = element_text(hjust = 0.5))


wilcox.test(Bifidobacterium_asteroides_noCOG_genes$mean[which(Bifidobacterium_asteroides_noCOG_genes$RVadj_agreement == "Non-sig.")],
            Bifidobacterium_asteroides_noCOG_genes$mean[which(Bifidobacterium_asteroides_noCOG_genes$RVadj_agreement == "Sig.")])



Gilliamella_apicola_noCOG_genes <- plasmid_prob_w_COG[which(plasmid_prob_w_COG$COG_category == "-" & plasmid_prob_w_COG$Species == "Gilliamella_apicola"), ]
Gilliamella_apicola_noCOG_genes$RVadj_agreement <- "Sig."

Gilliamella_apicola_noCOG_genes[which(rownames(Gilliamella_apicola_noCOG_genes) %in% putatively_high_mobilility_genes), "RVadj_agreement"] <- "Non-sig."

ggplot(data = Gilliamella_apicola_noCOG_genes, aes(x = RVadj_agreement, y = mean)) +
  geom_quasirandom(groupOnX = TRUE) +
  theme_bw() +
  ylab("Mean plasmid probability") +
  xlab("RVadj agreement with strain abundances") +
  ggtitle("Gilliamella apicola non-COG-associated genes") +
  theme(plot.title = element_text(hjust = 0.5))


wilcox.test(Gilliamella_apicola_noCOG_genes$mean[which(Gilliamella_apicola_noCOG_genes$RVadj_agreement == "Non-sig.")],
            Gilliamella_apicola_noCOG_genes$mean[which(Gilliamella_apicola_noCOG_genes$RVadj_agreement == "Sig.")])

