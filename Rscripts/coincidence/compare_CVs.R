### Compare coefficient of variation of genes in the same clusters for clusters thought to represent strain-specific groupings vs 
### co-occurence found only in breadth data (which could also represent primarily strain-specific variation).

rm(list = ls(all.names = TRUE))

compute_gene_cluster_CV <- function(gene_set, depth_table, presence_table, presence_cutoff, cluster_name) {
  # Compute CV of depth of all genes in a cluster based for all samples
  # where at least the specified proportion of the presence cutoff is met.
  
  presence_prop <- colSums(presence_table[gene_set, ]) / length(gene_set)
  
  present_samples <- names(presence_prop)[which(presence_prop >= presence_cutoff)]
  
  out_df <- data.frame(matrix(NA, nrow = length(present_samples), ncol = 5))
  colnames(out_df) <- c("cluster", "sample", "mean", "sd", "cv")
  
  out_df$cluster <- cluster_name
  out_df$sample <- present_samples
  
  if (length(present_samples) > 1) {
    out_df$mean <- sapply(depth_table[gene_set, present_samples], mean)
    out_df$sd <- sapply(depth_table[gene_set, present_samples], sd)
  } else {
    out_df$mean <- mean(depth_table[gene_set, present_samples])
    out_df$sd <- sd(depth_table[gene_set, present_samples])
  }
  
  out_df$cv <- out_df$sd / out_df$mean
  
  return(out_df)
  
}

Gilliamella_and_Snodgrassella_present_noncore <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_and_Snodgrassella_present_noncore.rds")

presence_cutoff <- 0.5
Gilliamella_and_Snodgrassella_present_noncore[Gilliamella_and_Snodgrassella_present_noncore >= presence_cutoff] <- 1
Gilliamella_and_Snodgrassella_present_noncore[Gilliamella_and_Snodgrassella_present_noncore < presence_cutoff] <- 0

Gilliamella_and_Snodgrassella_present_noncore <- Gilliamella_and_Snodgrassella_present_noncore[-which(rowSums(Gilliamella_and_Snodgrassella_present_noncore) == 0), ]
Gilliamella_and_Snodgrassella_present_noncore <- Gilliamella_and_Snodgrassella_present_noncore[-which(rowSums(Gilliamella_and_Snodgrassella_present_noncore) == 74), ]


Gilliamella_mean_depth <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_mean_depth_per_site.rds")
Snodgrassella_mean_depth <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Snodgrassella_mean_depth_per_site.rds")
mean_depth <- rbind(Gilliamella_mean_depth, Snodgrassella_mean_depth)
mean_depth <- mean_depth[-which(rowSums(mean_depth) == 0), ]

# Compute CV of depth for (likely) strain-level gene clusters.
Gilliamella_noncore_same_strain_cooccur_clusters <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Gilliamella_noncore_same_strain_cooccur_clusters.rds")
Snodgrassella_noncore_same_strain_cooccur_clusters <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Snodgrassella_noncore_same_strain_cooccur_clusters.rds")

names(Gilliamella_noncore_same_strain_cooccur_clusters) <- paste("Gilliamella_strain_cooccur", as.character(1:length(Gilliamella_noncore_same_strain_cooccur_clusters)), sep = "_")
names(Snodgrassella_noncore_same_strain_cooccur_clusters) <- paste("Snodgrassella_strain_cooccur", as.character(1:length(Snodgrassella_noncore_same_strain_cooccur_clusters)), sep = "_")

Gilliamella_strain_cococcur_CV <- lapply(X = names(Gilliamella_noncore_same_strain_cooccur_clusters),
              FUN = function(x) { 
          
                compute_gene_cluster_CV(gene_set = Gilliamella_noncore_same_strain_cooccur_clusters[[x]],
                                        depth_table = mean_depth,
                                        presence_table = Gilliamella_and_Snodgrassella_present_noncore,
                                        presence_cutoff = 0.75,
                                        cluster_name = x)
            
              })
Gilliamella_strain_cococcur_CV <- do.call(rbind, Gilliamella_strain_cococcur_CV)

Snodgrassella_strain_cococcur_CV <- lapply(X = names(Snodgrassella_noncore_same_strain_cooccur_clusters),
                                         FUN = function(x) { 
                                           
                                           compute_gene_cluster_CV(gene_set = Snodgrassella_noncore_same_strain_cooccur_clusters[[x]],
                                                                   depth_table = mean_depth,
                                                                   presence_table = Gilliamella_and_Snodgrassella_present_noncore,
                                                                   presence_cutoff = 0.75,
                                                                   cluster_name = x)
                                           
                                         })

Snodgrassella_strain_cococcur_CV <- do.call(rbind, Snodgrassella_strain_cococcur_CV)


# Now for other genes.
breadth_cooccur_clusters <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_cooccur_results_clustered_DBH0.3_components.rds")

breadth_cooccur_cluster_members <- list()
for (cluster_i in 1:breadth_cooccur_clusters$no) {
  breadth_cooccur_cluster_members[[cluster_i]] <- names(breadth_cooccur_clusters$membership)[which(breadth_cooccur_clusters$membership == cluster_i)]
}
names(breadth_cooccur_cluster_members) <- paste("breadth_cooccur", as.character(1:breadth_cooccur_clusters$no), sep = "_")


breadth_cococcur_CV <- lapply(X = names(breadth_cooccur_cluster_members),
                                         FUN = function(x) { 
                                           
                                           compute_gene_cluster_CV(gene_set = breadth_cooccur_cluster_members[[x]],
                                                                   depth_table = mean_depth,
                                                                   presence_table = Gilliamella_and_Snodgrassella_present_noncore,
                                                                   presence_cutoff = 0.75,
                                                                   cluster_name = x)
                                           
                                         })
breadth_cococcur_CV <- do.call(rbind, breadth_cococcur_CV)

par(mfrow = c(3, 1))
boxplot(cv ~ cluster, data = breadth_cococcur_CV, ylab = "CV in depth", main = "Breadth co-occurence clusters")
boxplot(cv ~ cluster, data = Gilliamella_strain_cococcur_CV, ylab = "CV in depth", main = "Gilliamella strain-content cooccurrence")
boxplot(cv ~ cluster, data = Snodgrassella_strain_cococcur_CV, ylab = "CV in depth", main = "Snodgrassella strain-content cooccurrence")
par(mfrow = c(1, 1))


boxplot(breadth_cococcur_CV$cv, Gilliamella_strain_cococcur_CV$cv, Snodgrassella_strain_cococcur_CV$cv,
        names = c("Breadth co-occurence clusters", "Gilliamella strain-content cooccurrence", "Snodgrassella strain-content cooccurrence"),
        ylab = "CV in depth")
