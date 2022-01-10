### After identifying gene pairs that cooccur more often in reference genomes and based on breadth profiles across MGS samples,
### figure out what clusters of genes there are that cooccur, which can be collapsed into more tractable profiles.

rm(list = ls(all.names = TRUE))

identify_coccuring_gene_clusters <- function(breadth_df) {
  
  gene_clusters = list()
  
  for (i in 1:nrow(breadth_df)) {
    
    if (breadth_df[i, "p.value"] >= 0.05) { next }
    
    gene1 <- breadth_df[i, "gene1"]
    gene2 <- breadth_df[i, "gene2"]
    
    cluster_matches <- c()
    
    if(length(gene_clusters) == 0) {
      gene_clusters[[1]] <- c(gene1, gene2)
      next
    }
    
    for (cluster_index in 1:length(gene_clusters)) {
      
      if ((gene1 %in% gene_clusters[[cluster_index]]) | (gene2 %in% gene_clusters[[cluster_index]])) {
        cluster_matches <- c(cluster_matches, cluster_index)
      }
    }
    
    if (length(cluster_matches) == 0) {
      # Define new cluster for these genes
      gene_clusters[[length(gene_clusters) + 1]] <- c(gene1, gene2)
      
    } else if (length(cluster_matches) == 1) {
      # Add these genes to an existing cluster
      gene_clusters[[cluster_matches[1]]] <- c(gene_clusters[[cluster_matches[1]]], gene1, gene2)
      gene_clusters[[cluster_matches[1]]] <- gene_clusters[[cluster_matches[1]]][-which(duplicated(gene_clusters[[cluster_matches[1]]]))]
      
    } else if (length(cluster_matches) > 1) {
      
      # Combine clusters that are linked by at least one of these genes.
      new_cluster <- c(gene1, gene2)
      
      for (j in cluster_matches) {
        new_cluster <- c(new_cluster, gene_clusters[[j]])
      }
      
      new_cluster <- new_cluster[-which(duplicated(new_cluster))]
      
      gene_clusters <- gene_clusters[-cluster_matches]
      
      gene_clusters[[length(gene_clusters) + 1]] <- new_cluster
      
    }
  }
  
  return(gene_clusters)
  
}

Gilliamella_noncore_same_strain_cooccur <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Gilliamella_noncore_same_strain_cooccur.rds")
Gilliamella_noncore_same_strain_cooccur_clusters <- identify_coccuring_gene_clusters(breadth_df = Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results)
saveRDS(object = Gilliamella_noncore_same_strain_cooccur_clusters, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Gilliamella_noncore_same_strain_cooccur_clusters.rds")

Snodgrassella_noncore_same_strain_cooccur <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Snodgrassella_noncore_same_strain_cooccur.rds")
Snodgrassella_noncore_same_strain_cooccur_clusters <- identify_coccuring_gene_clusters(breadth_df = Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results)
saveRDS(object = Snodgrassella_noncore_same_strain_cooccur_clusters, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Snodgrassella_noncore_same_strain_cooccur_clusters.rds")

