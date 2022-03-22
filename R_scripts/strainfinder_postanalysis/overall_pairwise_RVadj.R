# Compute RVadj between pairwise matrices: species, strain, and gene presence/absence
# (leaving out gene haplotypes of this analysis)

rm(list = ls(all.names = TRUE))

library(MatrixCorrelation)

species_presence <- readRDS(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/species_presence_50core_or_10percent.rds")

strain_abun <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/strain_abun.rds")

presence_matrix <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.matrix",
                              header = TRUE, sep = "\t", row.names = 1)

all_samples <- colnames(species_presence)
all_species <- rownames(species_presence)
presence_matrix <- t(presence_matrix[, all_samples])

category_data <- list()
category_data[["Species"]] <- t(species_presence)

for (species in names(strain_abun)) {
  
  species_strain_abun <- strain_abun[[species]]
  
  # Add in missing samples.
  missing_samples <- all_samples[which(! all_samples %in% rownames(species_strain_abun))]
  
  if (length(missing_samples) > 0) {
    missing_samples_df <- data.frame(matrix(0, nrow = length(missing_samples), ncol = ncol(species_strain_abun)))
    rownames(missing_samples_df) <- missing_samples
    colnames(missing_samples_df) <- colnames(species_strain_abun)
    species_strain_abun <- rbind(species_strain_abun, missing_samples_df)
  }

  category_data[[species]] <- as.matrix(species_strain_abun[all_samples, ])
  
  category_data[[paste("Genes:", species)]] <- presence_matrix[all_samples, grep(species, colnames(presence_matrix))]
  category_data[[paste("Genes:", species)]] <- as.matrix(category_data[[paste("Genes:", species)]][, which(colSums(category_data[[paste("Genes:", species)]]) > 0)])
  
}

all_categories <- c("Species", all_species, paste("Genes:", all_species, sep = " "))

RVadj_summary <- data.frame(matrix(NA, nrow = length(all_categories), ncol = length(all_categories)))
rownames(RVadj_summary) <- all_categories
colnames(RVadj_summary) <- all_categories

for (c1 in all_categories) {
 for (c2 in all_categories) {
   RVadj_summary[c1, c2] <- MatrixCorrelation::RVadj(category_data[[c1]],
                                                     category_data[[c2]])
 }
}



# Quick peek into this data.

library(ComplexHeatmap)

cor_break_cols = colorRamp2(c(0, 0.5, 0.95, 1), c("blue3", "white", "red3", "black"))

# Heatmap
heatmap <- draw(Heatmap(as.matrix(RVadj_summary), rect_gp = gpar(type = "none"), column_dend_side = "bottom",
                        name = "RVadj", row_names_side = "left",
                        cell_fun = function(j, i, x, y, w, h, fill) {
                          if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
                            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                          }
                        }))

heatmap
















# Do the same thing, but restrict to intersecting samples only for each case, as a comparison.

rm(list = ls(all.names = TRUE))

library(MatrixCorrelation)

species_presence <- readRDS(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/species_presence_50core_or_10percent.rds")

strain_abun <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/strain_abun.rds")

presence_matrix <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.matrix",
                              header = TRUE, sep = "\t", row.names = 1)

all_samples <- colnames(species_presence)
all_species <- rownames(species_presence)
presence_matrix <- t(presence_matrix[, all_samples])

category_data <- list()
category_data[["Species"]] <- t(species_presence)

for (species in names(strain_abun)) {
  
  species_strain_abun <- strain_abun[[species]]
  
  # Add in missing samples.
  missing_samples <- all_samples[which(! all_samples %in% rownames(species_strain_abun))]
  
  if (length(missing_samples) > 0) {
    missing_samples_df <- data.frame(matrix(0, nrow = length(missing_samples), ncol = ncol(species_strain_abun)))
    rownames(missing_samples_df) <- missing_samples
    colnames(missing_samples_df) <- colnames(species_strain_abun)
    species_strain_abun <- rbind(species_strain_abun, missing_samples_df)
  }
  
  category_data[[species]] <- as.matrix(species_strain_abun[all_samples, ])
  
  category_data[[paste("Genes:", species)]] <- presence_matrix[all_samples, grep(species, colnames(presence_matrix))]
  category_data[[paste("Genes:", species)]] <- as.matrix(category_data[[paste("Genes:", species)]][, which(colSums(category_data[[paste("Genes:", species)]]) > 0)])
  
}

all_categories <- c("Species", all_species, paste("Genes:", all_species, sep = " "))

RVadj_summary <- data.frame(matrix(NA, nrow = length(all_categories), ncol = length(all_categories)))
rownames(RVadj_summary) <- all_categories
colnames(RVadj_summary) <- all_categories

for (c1 in all_categories) {
  c1_samples <- rownames(category_data[[c1]])[which(rowSums(category_data[[c1]]) > 0)]
  
  for (c2 in all_categories) {
    c2_samples <- rownames(category_data[[c2]])[which(rowSums(category_data[[c2]]) > 0)]
    
    intersecting_samples <- c1_samples[which(c1_samples %in% c2_samples)]
    
    if (length(intersecting_samples) > 2) {
      RVadj_summary[c1, c2] <- MatrixCorrelation::RVadj(category_data[[c1]],
                                                        category_data[[c2]])
    }
  }
}



# Quick peek into this data.

library(ComplexHeatmap)


# Heatmap
heatmap <- draw(Heatmap(as.matrix(RVadj_summary), rect_gp = gpar(type = "none"), column_dend_side = "bottom",
                        name = "RVadj", row_names_side = "left",
                        cell_fun = function(j, i, x, y, w, h, fill) {
                          if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
                            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                          }
                        }))

heatmap


