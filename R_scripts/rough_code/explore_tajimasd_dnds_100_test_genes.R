rm(list = ls(all.names = TRUE))

phylotypes <- c("Bifidobacterium", "Firm4", "Firm5", "Gilliamella", "Snodgrassella")

tajimas_d_and_metric_files <- list.files(path = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/tajimas_d_and_metrics/",
                                         pattern = ".tsv",
                                         full.names = TRUE)

tajimas_d_and_metric_basenames <- list.files(path = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/tajimas_d_and_metrics/",
                                             pattern = ".tsv")

tajimas_d_and_metric_basenames <- gsub("_tajimas_d_and_metrics.tsv", "", tajimas_d_and_metric_basenames)

tajimas_d_and_metrics <- list()

for (i in 1:length(tajimas_d_and_metric_files)) {
  tajimas_d_and_metrics[[tajimas_d_and_metric_basenames[i]]] <- read.table(file = tajimas_d_and_metric_files[i], header = TRUE, stringsAsFactors = FALSE, sep = "\t")
}

metrics_summary <- list()

for (phylotype in phylotypes) {
  phylotype_test_genes <- grep(phylotype, tajimas_d_and_metric_basenames, value = TRUE)
  
  phylotype_metrics_df <- data.frame(matrix(NA, nrow = length(phylotype_test_genes), ncol = 4))
  rownames(phylotype_metrics_df) <- phylotype_test_genes
  colnames(phylotype_metrics_df) <- c("mean_D", "p_D_diff0", "mean_D_vs_all", "p_D_diff1_vs_all")
  
  for (gene in phylotype_test_genes) {
    gene_tajimas_d_and_metrics <- tajimas_d_and_metrics[[gene]]
    
    all_haplotype_D <- gene_tajimas_d_and_metrics[which(gene_tajimas_d_and_metrics$category == "all_haplotypes"), "D"]
      
    sample_D <- gene_tajimas_d_and_metrics[-which(gene_tajimas_d_and_metrics$category %in% c("all_haplotypes", "reference_seqs")), "D"]
    
    if (length(which(is.na(sample_D))) > 0) {
      sample_D <- sample_D[-which(is.na(sample_D))] 
    }
    
    if (length(sample_D) > 0) {
    
      phylotype_metrics_df[gene, c("mean_D", "p_D_diff0")] <- c(mean(sample_D), wilcox.test(sample_D, exact = FALSE)$p.value)
      
      if (! is.na(all_haplotype_D)) {
        
        sample_D_vs_all <- sample_D / all_haplotype_D
        
        phylotype_metrics_df[gene, c("mean_D_vs_all", "p_D_diff1_vs_all")] <- c(mean(sample_D_vs_all), wilcox.test(sample_D_vs_all, mu = 1, exact = FALSE)$p.value) 
      }
      
    }
  }
  
  metrics_summary[[phylotype]] <- phylotype_metrics_df
}


par(mfrow=c(1, 2))

boxplot(metrics_summary$Bifidobacterium$mean_D,
        metrics_summary$Firm4$mean_D,
        metrics_summary$Firm5$mean_D,
        metrics_summary$Gilliamella$mean_D,
        metrics_summary$Snodgrassella$mean_D,
        names = phylotypes,
        main = "Mean sample D")


boxplot(metrics_summary$Bifidobacterium$mean_D_vs_all,
        metrics_summary$Firm4$mean_D_vs_all,
        metrics_summary$Firm5$mean_D_vs_all,
        metrics_summary$Gilliamella$mean_D_vs_all,
        metrics_summary$Snodgrassella$mean_D_vs_all,
        names = phylotypes,
        main = "Mean of sample D / all haplotype D")











dnds_files <- list.files(path = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/dnds_out/",
                                         pattern = "_dnds.tsv",
                                         full.names = TRUE)

dnds_files <- dnds_files[-grep("haplotype_dnds.tsv", dnds_files)]

dnds_basenames <- list.files(path = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/dnds_out/",
                         pattern = "_dnds.tsv",
                         full.names = FALSE)

dnds_basenames <- dnds_basenames[-grep("haplotype_dnds.tsv", dnds_basenames)]

dnds_basenames <- gsub("_dnds.tsv", "", dnds_basenames)





dnds <- list()

for (i in 1:length(dnds_files)) {
  dnds[[dnds_basenames[i]]] <- read.table(file = dnds_files[i], header = TRUE, stringsAsFactors = FALSE, sep = "\t")
}

dnds_summary <- list()

for (phylotype in phylotypes) {
  phylotype_test_genes <- grep(phylotype, dnds_basenames, value = TRUE)
  
  phylotype_dnds_df <- data.frame(matrix(NA, nrow = length(phylotype_test_genes), ncol = 4))
  rownames(phylotype_dnds_df) <- phylotype_test_genes
  colnames(phylotype_dnds_df) <- c("mean_dnds", "p_dnds_diff0", "mean_dnds_vs_all", "p_dnds_diff1_vs_all")
  
  for (gene in phylotype_test_genes) {
    gene_dnds <- dnds[[gene]]
    
    all_haplotype_dnds <- gene_dnds[which(gene_dnds$category == "all_haplotypes"), "dnds"]
    
    sample_dnds <- gene_dnds[-which(gene_dnds$category %in% c("all_haplotypes", "reference_seqs")), "dnds"]
    
    if (length(which(is.na(sample_dnds))) > 0) {
      sample_dnds <- sample_dnds[-which(is.na(sample_dnds))] 
    }
    
    if (length(sample_dnds) > 0) {
      
      phylotype_dnds_df[gene, c("mean_dnds", "p_dnds_diff0")] <- c(mean(sample_dnds), wilcox.test(sample_dnds, exact = FALSE)$p.value)
      
      if (! is.na(all_haplotype_dnds)) {
        
        sample_dnds_vs_all <- sample_dnds - all_haplotype_dnds
        
        phylotype_dnds_df[gene, c("mean_dnds_vs_all", "p_dnds_diff1_vs_all")] <- c(mean(sample_dnds_vs_all), wilcox.test(sample_dnds_vs_all, mu = 1, exact = FALSE)$p.value) 
      }
      
    }
  }
  
  dnds_summary[[phylotype]] <- phylotype_dnds_df
}


boxplot(dnds_summary$Bifidobacterium$mean_dnds,
        dnds_summary$Firm4$mean_dnds,
        dnds_summary$Firm5$mean_dnds,
        dnds_summary$Gilliamella$mean_dnds,
        dnds_summary$Snodgrassella$mean_dnds,
        names = phylotypes,
        main = "Mean dN/dS")


boxplot(dnds_summary$Bifidobacterium$mean_dnds_vs_all,
        dnds_summary$Firm4$mean_dnds_vs_all,
        dnds_summary$Firm5$mean_dnds_vs_all,
        dnds_summary$Gilliamella$mean_dnds_vs_all,
        dnds_summary$Snodgrassella$mean_dnds_vs_all,
        names = phylotypes,
        main = "Mean dN/dS across all samples - dN/dS of all haplotypes")

