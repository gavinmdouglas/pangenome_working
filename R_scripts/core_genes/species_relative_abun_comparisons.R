### Compare relative abundances at species-level based on mean depth of genes with variants.
### Mainly motivated by the question of how the relative abundances Gilliamella apicola and Bartonella apis relate
### (i.e., in addition to their presence/absence profiles showing a negative association, is the same true for their)
### relative abundance?

### Note that I looked at rel abundance based on cases where sp was called as present only and also for cases when it was below cut-off (species_relabun_all_genes)

rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(kableExtra)
library(reshape2)
library(ggplot2)

pandora_output <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.matrix",
                             header = TRUE, sep = "\t", row.names = 1)

variant_mean_depth <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_mean_variant_site_depth.tsv",
                                 header = TRUE, sep = "\t", row.names = 1)
variant_mean_depth[variant_mean_depth == 0] <- NA

species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species_present.txt", stringsAsFactors = FALSE)$V1

species_presence <- data.frame(matrix(0, nrow = length(species), ncol = ncol(pandora_output)))
colnames(species_presence) <- colnames(pandora_output)
rownames(species_presence) <- species

species_relabun <- data.frame(matrix(0, nrow = length(species), ncol = ncol(pandora_output)))
colnames(species_relabun) <- colnames(pandora_output)
rownames(species_relabun) <- species

species_relabun_all_genes <- data.frame(matrix(0, nrow = length(species), ncol = ncol(pandora_output)))
colnames(species_relabun_all_genes) <- colnames(pandora_output)
rownames(species_relabun_all_genes) <- species

panaroo_only_potential_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/panaroo_only_potential_core.rds")
checkm_panaroo_potential_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/checkm_panaroo_passed_potential_core.rds")

core_genes <- c(panaroo_only_potential_core, checkm_panaroo_potential_core)
core_genes <- core_genes[species]
#saveRDS(object = core_genes,
#        file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/species_present_all_potential_core.rds")

potential_core_present <- list()

any_present <- c()

for (sp in species) {
  potential_core_present[[sp]] <- core_genes[[sp]][which(core_genes[[sp]] %in% rownames(pandora_output))]
  any_present <- c(any_present, length(grep(sp, rownames(pandora_output))))
  
  num_core_present <- colSums(pandora_output[potential_core_present[[sp]], ])
  
  samples_i_with_sp <- which(num_core_present >= 10 | (num_core_present / length(core_genes[[sp]])) >= 0.1)
  
  if (length(samples_i_with_sp) > 0) {
    species_presence[sp, samples_i_with_sp] <- 1
  }
  
  core_genes_w_depth <- potential_core_present[[sp]][which(potential_core_present[[sp]] %in% rownames(variant_mean_depth))]
  species_relabun_all_genes[sp, ] <- colMeans(variant_mean_depth[core_genes_w_depth, colnames(species_relabun_all_genes)], na.rm = TRUE)
  species_relabun[sp, samples_i_with_sp] <- colMeans(variant_mean_depth[core_genes_w_depth, colnames(species_relabun)[samples_i_with_sp]], na.rm = TRUE)
}

species_relabun[is.na(species_relabun)] <- 0

species_relabun <- data.frame(sweep(species_relabun, 2, colSums(species_relabun), "/")) * 100

species_relabun_t <- data.frame(t(species_relabun))

plot(species_relabun_t$Gilliamella_apicola, species_relabun_t$Bartonella_apis,
     xlab = "Gilliamella apicola rel. abun. (%)", ylab = "Bartonella apis rel. abun. (%)")

Heatmap(as.matrix(species_presence),
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        clustering_distance_rows = "binary",
        clustering_distance_columns = "binary")


# Stacked barchart.

quick_hclust <- hclust(dist(t(species_relabun)))
quick_hclust$labels[quick_hclust$order]

species_relabun_long <- melt(as.matrix(species_relabun))
species_relabun_long <- species_relabun_long[-which(species_relabun_long$value == 0), ]

colnames(species_relabun_long)[1] <- "Species"

species_relabun_long$Var2 <- factor(species_relabun_long$Var2, levels = quick_hclust$labels[quick_hclust$order])

plot_colours <- c("#58c07b",
                   "#593789",
                   "#91b23e",
                   "#c170c8",
                   "#56903e",
                   "#6d80d8",
                   "#cd9c2e",
                   "#b6518f",
                   "#43c9b0",
                   "#d6565a",
                   "#9f9047",
                   "#bb486a",
                   "#ba9940",
                   "#973f2b",
                   "#c86938")

ggplot(data = species_relabun_long, aes(y = Var2, x = value, fill = Species)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = plot_colours) +
  xlab("Relative abundance (%)") +
  ylab("Sample") +
  theme_bw()
