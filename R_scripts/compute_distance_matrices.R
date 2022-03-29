# Compute distance matrices based on species, strain, and gene presence / absence.

# Compute presence/absence distance based on Jaccard and "alpha" approach.

# Also produce distance matrices for strains based on Bray-Curtis. 

rm(list = ls(all.names = TRUE))

library(vegan)
library(MatrixCorrelation)


ML.Alpha_over_matrix <- function(in_matrix) {
  
  in_matrix <- in_matrix[, which(colSums(in_matrix) > 0)]
  
  if (is.null(rownames(in_matrix))) { stop("Rownames needed.") }
  
  matrix_samples <- rownames(in_matrix)

  alpha_output <- data.frame(matrix(0, nrow = nrow(in_matrix), ncol = nrow(in_matrix)))
  rownames(alpha_output) <- matrix_samples
  colnames(alpha_output) <- matrix_samples
  
  num_features <- ncol(in_matrix)
  
  for (s1 in matrix_samples) {
    
    s1_i <- which(matrix_samples == s1)
    
    s1_present <- which(as.numeric(in_matrix[s1, ]) > 0)
    
    num_s1_present <- length(s1_present)
    
    if (num_s1_present == 0) { next }
    
    # Include comparison with self as a sanity check for now.
    remaining_samples <- matrix_samples[s1_i:length(matrix_samples)]
    
    for (s2 in remaining_samples) {

      s2_present <- which(as.numeric(in_matrix[s2, ]) > 0)
      
      intersecting <- length(which(s1_present %in% s2_present))
      
      num_s2_present <- length(s2_present)
      
      if (num_s2_present == 0) { next }
      
      alpha_output[s1, s2] <- ML.Alpha(intersecting, c(num_s1_present, num_s2_present, num_features))$est
      
    }
  }

  return(alpha_output)
}

species_presence <- readRDS(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/species_presence_50core_or_10percent.rds")

strain_abun <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/strain_abun.rds")

presence_matrix <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.matrix",
                              header = TRUE, sep = "\t", row.names = 1)


all_samples <- colnames(species_presence)

presence_matrix <- presence_matrix[, all_samples]

jaccard_dist_matrices <- list()

jaccard_dist_matrices[["Species"]] <- as.matrix(dist(x = t(species_presence), method = "binary", diag = TRUE))

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
  
  species_strain_abun <- species_strain_abun[all_samples, ]

  jaccard_dist_matrices[[species]] <- as.matrix(dist(x = species_strain_abun, method = "binary", diag = TRUE))
  
  presence_matrix_species <- presence_matrix[grep(species, rownames(presence_matrix)), ]
  jaccard_dist_matrices[[paste(species, "genes")]] <- as.matrix(dist(x = t(presence_matrix_species), method = "binary", diag = TRUE))
  
}

tmp1 <- ML.Alpha_over_matrix(t(presence_matrix_species))
tmp2 <- ML.Alpha_over_matrix(species_strain_abun)

mantel(tmp1, tmp2, method="pearson", permutations=999, strata = NULL,
       na.rm = FALSE, parallel = getOption("mc.cores"))

species_strain_abun_rand <- species_strain_abun[sample(1:nrow(species_strain_abun)), ]

allCorrelations(as.matrix(t(presence_matrix_species)),
                as.matrix(species_strain_abun),
                methods = "RVadj", ncomp1 = 5, ncomp2 = 5)

MatrixCorrelation::SMI(as.matrix(t(presence_matrix_species)),
                         as.matrix(species_strain_abun))

tmp1_subset <- tmp1[! is.na(tmp1)]
tmp2_subset <- tmp2[! is.na(tmp2)]

mantel(jaccard_dist_matrices[[paste(species, "genes")]],
       jaccard_dist_matrices[[species]],
       method="pearson", permutations=999, strata = NULL,
       na.rm = FALSE, parallel = getOption("mc.cores"))

jaccard_dist_matrices


