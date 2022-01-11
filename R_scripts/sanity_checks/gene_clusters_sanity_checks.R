### Run sanity checks on gene clusters that were produced.

rm(list = ls(all.names = TRUE))

Gilliamella_noncore_same_strain_cooccur <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Gilliamella_noncore_same_strain_cooccur.rds")
Gilliamella_noncore_same_strain_cooccur_clusters <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Gilliamella_noncore_same_strain_cooccur_clusters.rds")

Snodgrassella_noncore_same_strain_cooccur <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Snodgrassella_noncore_same_strain_cooccur.rds")
Snodgrassella_noncore_same_strain_cooccur_clusters <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Snodgrassella_noncore_same_strain_cooccur_clusters.rds")


# Check that the number of unique genes in table and cluster list matches.
unique_Gilliamella_genes <- c(Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results$gene1, Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results$gene2)
unique_Gilliamella_genes <- unique_Gilliamella_genes[-which(duplicated(unique_Gilliamella_genes))]
length(unique_Gilliamella_genes) == sum(sapply(Gilliamella_noncore_same_strain_cooccur_clusters, length))

unique_Snodgrassella_genes <- c(Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results$gene1, Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results$gene2)
unique_Snodgrassella_genes <- unique_Snodgrassella_genes[-which(duplicated(unique_Snodgrassella_genes))]
length(unique_Snodgrassella_genes) == sum(sapply(Snodgrassella_noncore_same_strain_cooccur_clusters, length))



# Check for a couple of frequent genes that all pairing genes are in the same cluster.

identify_matching_cluster <- function(gene, cluster_list) {

  for (i in 1:length(cluster_list)) {
    if (gene %in% cluster_list[[i]]) {
      return(i) 
    }
  }
  return(NA)
}

# Gilli_group_2155
Gilli_group_2155_matches <- c(Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results[which(Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results$gene1 == "Gilli_group_2155"), "gene2"],
                              Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results[which(Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results$gene2 == "Gilli_group_2155"), "gene1"])
table(sapply(X = c("Gilli_group_2155", Gilli_group_2155_matches), FUN = identify_matching_cluster, cluster_list = Gilliamella_noncore_same_strain_cooccur_clusters))


# Gilli_group_2137
Gilli_group_2137_matches <- c(Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results[which(Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results$gene1 == "Gilli_group_2137"), "gene2"],
                              Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results[which(Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results$gene2 == "Gilli_group_2137"), "gene1"])
table(sapply(X = c("Gilli_group_2137", Gilli_group_2137_matches), FUN = identify_matching_cluster, cluster_list = Gilliamella_noncore_same_strain_cooccur_clusters))


# Snod_ogt
Snod_ogt_matches <- c(Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results[which(Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results$gene1 == "Snod_ogt"), "gene2"],
                              Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results[which(Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results$gene2 == "Snod_ogt"), "gene1"])
table(sapply(X = c("Snod_ogt", Snod_ogt_matches), FUN = identify_matching_cluster, cluster_list = Snodgrassella_noncore_same_strain_cooccur_clusters))


# Snod_group_43
Snod_group_43_matches <- c(Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results[which(Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results$gene1 == "Snod_group_43"), "gene2"],
                      Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results[which(Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results$gene2 == "Snod_group_43"), "gene1"])
table(sapply(X = c("Snod_group_43", Snod_group_43_matches), FUN = identify_matching_cluster, cluster_list = Snodgrassella_noncore_same_strain_cooccur_clusters))



# Check for a couple of clusters with very few genes that they are correct.

Gilli_cluster1_genes <- Gilliamella_noncore_same_strain_cooccur_clusters[[1]]
Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results[which(Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results$gene1 %in% Gilli_cluster1_genes | Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results$gene2 %in% Gilli_cluster1_genes), ]

Gilli_cluster2_genes <- Gilliamella_noncore_same_strain_cooccur_clusters[[2]]
Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results[which(Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results$gene1 %in% Gilli_cluster2_genes | Gilliamella_noncore_same_strain_cooccur$breadth_coccur_results$gene2 %in% Gilli_cluster2_genes), ]


Snod_cluster1_genes <- Snodgrassella_noncore_same_strain_cooccur_clusters[[1]]
Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results[which(Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results$gene1 %in% Snod_cluster1_genes | Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results$gene2 %in% Snod_cluster1_genes), ]

Snod_cluster3_genes <- Snodgrassella_noncore_same_strain_cooccur_clusters[[3]]
Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results[which(Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results$gene1 %in% Snod_cluster3_genes | Snodgrassella_noncore_same_strain_cooccur$breadth_coccur_results$gene2 %in% Snod_cluster3_genes), ]

