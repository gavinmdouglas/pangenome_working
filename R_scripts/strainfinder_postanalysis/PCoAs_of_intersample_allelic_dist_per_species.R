rm(list = ls(all.names = TRUE))

library(ape)
library(vegan)
library(ggplot2)
library(ggbeeswarm)
library(reshape2)

# Compute PCs based on mean distance matrix of all other species.
combine_other_dists <- function(sp2ignore, metric, all_dists) {
  
  focal_matrix <- as.matrix(all_dists[[sp]][[m]])
  
  focal_samples <- rownames(focal_matrix)
  num_focal_samples <- length(focal_samples)
  
  summed_dists <- data.frame(matrix(0, nrow = num_focal_samples, ncol = num_focal_samples))
  rownames(summed_dists) <- focal_samples
  colnames(summed_dists) <- focal_samples
  
  num_comparisons <- data.frame(matrix(0, nrow = num_focal_samples, ncol = num_focal_samples))
  rownames(num_comparisons) <- focal_samples
  colnames(num_comparisons) <- focal_samples
  
  other_sp <- names(all_dists)[which(names(all_dists) != sp2ignore)]
  
  for (other in other_sp) {
    
    other_dist <- as.matrix(all_dists[[other]][[m]])
    intersecting_samples <- rownames(other_dist)[which(rownames(other_dist) %in% focal_samples)]
    
    if (length(intersecting_samples) > 0) {
      
      other_dist <- other_dist[intersecting_samples, intersecting_samples, drop = FALSE]
      
      other_dist_nonzero <- other_dist
      other_dist_nonzero[other_dist_nonzero > 0] <- 1
      
      summed_dists[intersecting_samples, intersecting_samples] <- summed_dists[intersecting_samples, intersecting_samples] + other_dist
      num_comparisons[intersecting_samples, intersecting_samples] <- num_comparisons[intersecting_samples, intersecting_samples] + other_dist_nonzero
      
    }
    
  }
  
  mean_other_dist <- summed_dists / num_comparisons
  mean_other_dist[is.na(mean_other_dist)] <- 0
  
  return(mean_other_dist)
  
}

sample_metadata <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard.2019.2020_metadata.tsv",
                              row.names = 1, header = TRUE, sep = "\t")

per_sp_dists <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/allelic_mean_intersample_dists.rds")

metrics <-  c("jaccard", "simpson", "jsd")

adonis_summary <- data.frame(matrix(NA, nrow = length(per_sp_dists) * 15, ncol = 11))
colnames(adonis_summary) <- c("Species", "Metric",
                              "Country_R2", "Country_P",
                              "Apiary_R2", "Apiary_P",
                              "Year_R2", "Year_P",
                              "Age_R2", "Age_P",
                              "Residuals_R2")

adonis_summary$Species <- rep(names(per_sp_dists$jaccard), each = 3)
adonis_summary$Metric <- rep(metrics, 15)

pcoa_plots <- list()

dist_matrices <- list()

for (sp in names(per_sp_dists$jaccard)) {
  
  pcoa_plots[[sp]] <- list()
  
  dist_matrices[[sp]] <- list()
  
  for (m in metrics) {
  
    dist_matrix <- per_sp_dists[[m]][[sp]]
    
    dist_matrix_rownames <- rownames(as.matrix(dist_matrix))
    
    sample_metadata_subset <- sample_metadata[dist_matrix_rownames, ]
    
    dist_matrices[[sp]][[m]] <- dist_matrix
    
    pcoa_out <- pcoa(dist_matrix, correction = "cailliez")
    
    if (ncol(pcoa_out$vectors) > 1) {

      axes_out <- cbind(data.frame(x1 = pcoa_out$vectors[, 1],
                                   x2 = pcoa_out$vectors[, 2]),
                        sample_metadata_subset)
      
      if ("Relative_eig" %in% colnames(pcoa_out$values)) {
      
        pcoa_plots[[sp]][[m]] <- ggplot(data = axes_out, aes(x = x1, y = x2, colour = Country, shape = Apiary)) +
                                        geom_point(size = 3) +
                                        theme_bw() +
                                        xlab(paste("PC1 (", as.character(round(pcoa_out$values$Relative_eig[1], 3)), "%)", sep = "")) +
                                        ylab(paste("PC2 (", as.character(round(pcoa_out$values$Relative_eig[2], 3)), "%)", sep = "")) +
                                        ggtitle(paste(sp, m, sep = " - ")) +
                                        theme(plot.title = element_text(hjust = 0.5))
        
      } else if ("Rel_corr_eig" %in% colnames(pcoa_out$values)) {
        
        pcoa_plots[[sp]][[m]] <- ggplot(data = axes_out, aes(x = x1, y = x2, colour = Country, shape = Apiary)) +
                                        geom_point(size = 3) +
                                        theme_bw() +
                                        xlab(paste("PC1 (", as.character(round(pcoa_out$values$Rel_corr_eig[1], 3)), "%)", sep = "")) +
                                        ylab(paste("PC2 (", as.character(round(pcoa_out$values$Rel_corr_eig[2], 3)), "%)", sep = "")) +
                                        ggtitle(paste(sp, m, sep = " - ")) +
                                        theme(plot.title = element_text(hjust = 0.5))
        
      } else {
       
        stop("Eigenvalue percent column not found.") 
        
      }
      
    }
    
    if (length(dist_matrix_rownames) >= 10) {
      
      # Only test factors for which there are at least 5 samples in at least two different groups.
      factors_to_include <- c()
      for (fac in c("Country", "Apiary", "Year", "Age")) {
        fac_breakdown <- table(sample_metadata_subset[, fac])
        if (length(which(fac_breakdown >= 5)) > 1) {
          factors_to_include <- c(factors_to_include, fac)
        }
      }
      
      if (length(factors_to_include) == 0) { next }
      
      adonis_formula <- as.formula(paste("dist_matrix ~ ", paste(factors_to_include, collapse = " + "), sep = ""))
      
      adonis_output <- data.frame(adonis(formula = adonis_formula, data = sample_metadata_subset, permutations = 9999)$aov.tab)
      
      row_i <- which(adonis_summary$Species == sp & adonis_summary$Metric == m)
      
      if ("Country" %in% rownames(adonis_output)) {
        adonis_summary[row_i, c("Country_R2", "Country_P")] <- as.numeric(adonis_output["Country", c("R2", "Pr..F.")])
      }
      
      if ("Apiary" %in% rownames(adonis_output)) {
        adonis_summary[row_i, c("Apiary_R2", "Apiary_P")] <- as.numeric(adonis_output["Apiary", c("R2", "Pr..F.")])
      }
      
      if ("Year" %in% rownames(adonis_output)) {
        adonis_summary[row_i, c("Year_R2", "Year_P")] <- as.numeric(adonis_output["Year", c("R2", "Pr..F.")])
      }
      
      if ("Age" %in% rownames(adonis_output)) {
        adonis_summary[row_i, c("Age_R2", "Age_P")] <- as.numeric(adonis_output["Age", c("R2", "Pr..F.")])
      }
 
      adonis_summary[row_i, "Residuals_R2"] <- adonis_output["Residuals", "R2"]
      
    }
    
  }
  
}

adonis_summary$Country_R2[which(adonis_summary$Country_P >= 0.05)] <- NA
adonis_summary$Apiary_R2[which(adonis_summary$Apiary_P >= 0.05)] <- NA
adonis_summary$Year_R2[which(adonis_summary$Year_P >= 0.05)] <- NA
adonis_summary$Age_R2[which(adonis_summary$Age_P >= 0.05)] <- NA

adonis_summary <- adonis_summary[, grep("_P$", colnames(adonis_summary), value = TRUE, invert = TRUE)]

adonis_summary_long <- melt(adonis_summary)

adonis_summary_long$variable <- gsub("_R2$", "", adonis_summary_long$variable)

adonis_summary_long <- adonis_summary_long[which(! is.na(adonis_summary_long$value)), ]

adonis_summary_long$Metric[which(adonis_summary_long$Metric == "jsd")] <- "Jensen-Shannon"
adonis_summary_long$Metric[which(adonis_summary_long$Metric == "jaccard")] <- "Jaccard"
adonis_summary_long$Metric[which(adonis_summary_long$Metric == "simpson")] <- "Simpson"

adonis_summary_long_no_residual <- adonis_summary_long[-which(adonis_summary_long$variable == "Residuals"), ]

adonis_summary_long_no_residual$variable <- factor(adonis_summary_long_no_residual$variable,
                                                   levels = c("Country", "Apiary", "Year", "Age"))

ggplot(data = adonis_summary_long_no_residual, aes(x = value, y = Species, colour = Metric)) +
  geom_point(size = 8, alpha = 0.5) +
  facet_wrap(. ~ variable) +
  scale_colour_manual(values = c("burlywood2", "brown", "chartreuse")) +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 1)) +
  xlab(expression(paste("PERMANOVA R"^"2"))) +
  ylab("")


adonis_summary_long_residual <- adonis_summary_long[which(adonis_summary_long$variable == "Residuals"), ]

ggplot(data = adonis_summary_long_residual, aes(x = value, y = Species, colour = Metric)) +
  geom_point(size = 8, alpha = 0.5) +
  scale_colour_manual(values = c("burlywood2", "brown", "chartreuse")) +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 1)) +
  xlab(expression(paste("PERMANOVA R"^"2"))) +
  ggtitle("Residuals") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("")















# Re-run PERMANOVAs but with first three PCs based on mean distance matrices of **all other species** included.

adonis_summary_w_PCs <- data.frame(matrix(NA, nrow = length(per_sp_dists) * 15, ncol = 17))
colnames(adonis_summary_w_PCs) <- c("Species", "Metric",
                                    "Country_R2", "Country_P",
                                    "Apiary_R2", "Apiary_P",
                                    "Year_R2", "Year_P",
                                    "Age_R2", "Age_P",
                                    "PC1_R2", "PC1_P",
                                    "PC2_R2", "PC2_P",
                                    "PC3_R2", "PC3_P",
                                    "Residuals_R2")

adonis_summary_w_PCs$Species <- rep(names(per_sp_dists$jaccard), each = 3)
adonis_summary_w_PCs$Metric <- rep(metrics, 15)

for (sp in names(per_sp_dists$jaccard)) {
  
  for (m in metrics) {
    
    dist_matrix <- dist_matrices[[sp]][[m]]
    
    sample_metadata_subset <- sample_metadata[rownames(as.matrix(dist_matrix)), ]
    
    pcoa_out <- pcoa(dist_matrix, correction = "cailliez")
    
    if (ncol(as.matrix(dist_matrix)) >= 10) {
      
      # Only test factors for which there are at least 5 samples in at least two different groups.
      factors_to_include <- c()
      for (fac in c("Country", "Apiary", "Year", "Age")) {
        fac_breakdown <- table(sample_metadata_subset[, fac])
        if (length(which(fac_breakdown >= 5)) > 1) {
          factors_to_include <- c(factors_to_include, fac)
        }
      }
      
      mean_other_dist <- combine_other_dists(sp2ignore=sp, metric=m,
                                             all_dists=dist_matrices)
      
      mean_other_pcoa <-  pcoa(dist_matrix, correction = "cailliez")
      
      if (ncol(mean_other_pcoa$vectors) > 1) {
        
        sample_metadata_subset <- cbind(data.frame(Other_PC1 = mean_other_pcoa$vectors[, 1],
                                                   Other_PC2 = mean_other_pcoa$vectors[, 2],
                                                   Other_PC3 = mean_other_pcoa$vectors[, 3]),
                                        sample_metadata_subset)
        
        factors_to_include <- c(factors_to_include, c("Other_PC1", "Other_PC2", "Other_PC3"))
      }
      
      
      if (length(factors_to_include) == 0) { next }
      
      adonis_formula <- as.formula(paste("dist_matrix ~ ", paste(factors_to_include, collapse = " + "), sep = ""))
      
      adonis_output <- data.frame(adonis(formula = adonis_formula, data = sample_metadata_subset, permutations = 9999)$aov.tab)
      
      row_i <- which(adonis_summary_w_PCs$Species == sp & adonis_summary_w_PCs$Metric == m)
      
      if ("Country" %in% rownames(adonis_output)) {
        adonis_summary_w_PCs[row_i, c("Country_R2", "Country_P")] <- as.numeric(adonis_output["Country", c("R2", "Pr..F.")])
      }
      
      if ("Apiary" %in% rownames(adonis_output)) {
        adonis_summary_w_PCs[row_i, c("Apiary_R2", "Apiary_P")] <- as.numeric(adonis_output["Apiary", c("R2", "Pr..F.")])
      }
      
      if ("Year" %in% rownames(adonis_output)) {
        adonis_summary_w_PCs[row_i, c("Year_R2", "Year_P")] <- as.numeric(adonis_output["Year", c("R2", "Pr..F.")])
      }
      
      if ("Age" %in% rownames(adonis_output)) {
        adonis_summary_w_PCs[row_i, c("Age_R2", "Age_P")] <- as.numeric(adonis_output["Age", c("R2", "Pr..F.")])
      }
      
      if ("Other_PC1" %in% rownames(adonis_output)) {
        adonis_summary_w_PCs[row_i, c("PC1_R2", "PC1_P")] <- as.numeric(adonis_output["Other_PC1", c("R2", "Pr..F.")])
        adonis_summary_w_PCs[row_i, c("PC2_R2", "PC2_P")] <- as.numeric(adonis_output["Other_PC2", c("R2", "Pr..F.")])
        adonis_summary_w_PCs[row_i, c("PC3_R2", "PC3_P")] <- as.numeric(adonis_output["Other_PC3", c("R2", "Pr..F.")])
      }
      
      adonis_summary_w_PCs[row_i, "Residuals_R2"] <- adonis_output["Residuals", "R2"]
      
    }
    
  }
  
}

adonis_summary_w_PCs$Country_R2[which(adonis_summary_w_PCs$Country_P >= 0.05)] <- NA
adonis_summary_w_PCs$Apiary_R2[which(adonis_summary_w_PCs$Apiary_P >= 0.05)] <- NA
adonis_summary_w_PCs$Year_R2[which(adonis_summary_w_PCs$Year_P >= 0.05)] <- NA
adonis_summary_w_PCs$Age_R2[which(adonis_summary_w_PCs$Age_P >= 0.05)] <- NA
adonis_summary_w_PCs$PC1_R2[which(adonis_summary_w_PCs$PC1_P >= 0.05)] <- NA
adonis_summary_w_PCs$PC2_R2[which(adonis_summary_w_PCs$PC2_P >= 0.05)] <- NA
adonis_summary_w_PCs$PC3_R2[which(adonis_summary_w_PCs$PC3_P >= 0.05)] <- NA

adonis_summary_w_PCs <- adonis_summary_w_PCs[, grep("_P$", colnames(adonis_summary_w_PCs), value = TRUE, invert = TRUE)]

adonis_summary_w_PCs_long <- melt(adonis_summary_w_PCs)

adonis_summary_w_PCs_long$variable <- gsub("_R2$", "", adonis_summary_w_PCs_long$variable)

adonis_summary_w_PCs_long <- adonis_summary_w_PCs_long[which(! is.na(adonis_summary_w_PCs_long$value)), ]

adonis_summary_w_PCs_long$Metric[which(adonis_summary_w_PCs_long$Metric == "jsd")] <- "Jensen-Shannon"
adonis_summary_w_PCs_long$Metric[which(adonis_summary_w_PCs_long$Metric == "jaccard")] <- "Jaccard"
adonis_summary_w_PCs_long$Metric[which(adonis_summary_w_PCs_long$Metric == "simpson")] <- "Simpson"

adonis_summary_w_PCs_long_no_residual <- adonis_summary_w_PCs_long[-which(adonis_summary_w_PCs_long$variable == "Residuals"), ]

adonis_summary_w_PCs_long_no_residual_PC_only <- adonis_summary_w_PCs_long_no_residual[grep("PC", adonis_summary_w_PCs_long_no_residual$variable), ]

adonis_summary_w_PCs_long_no_residual_PC_only$variable <- factor(adonis_summary_w_PCs_long_no_residual_PC_only$variable,
                                                   levels = c("PC1", "PC2", "PC3"))

ggplot(data = adonis_summary_w_PCs_long_no_residual_PC_only, aes(x = value, y = Species, colour = Metric)) +
  geom_point(size = 8, alpha = 0.5) +
  facet_wrap(. ~ variable) +
  scale_colour_manual(values = c("burlywood2", "brown", "chartreuse")) +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 1)) +
  xlab(expression(paste("PERMANOVA R"^"2"))) +
  ylab("")


adonis_summary_w_PCs_long_residual <- adonis_summary_w_PCs_long[which(adonis_summary_w_PCs_long$variable == "Residuals"), ]

adonis_summary_w_PCs_long_residual$value[which(adonis_summary_w_PCs_long_residual$value < 0)] <- 0

ggplot(data = adonis_summary_w_PCs_long_residual, aes(x = value, y = Species, colour = Metric)) +
  geom_point(size = 8, alpha = 0.5) +
  scale_colour_manual(values = c("burlywood2", "brown", "chartreuse")) +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 1)) +
  xlab(expression(paste("PERMANOVA R"^"2"))) +
  ggtitle("Residuals") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("")


