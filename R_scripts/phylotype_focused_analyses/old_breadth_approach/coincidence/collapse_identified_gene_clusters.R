### Make new breadth tables with gene clusters collapsed into single categories.

rm(list = ls(all.names = TRUE))

Gilliamella_and_Snodgrassella_present_noncore <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_and_Snodgrassella_present_noncore.rds")
Gilliamella_and_Snodgrassella_present_noncore <- Gilliamella_and_Snodgrassella_present_noncore[-which(rowSums(Gilliamella_and_Snodgrassella_present_noncore) == 0), ]

Gilliamella_noncore_same_strain_cooccur_clusters <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Gilliamella_noncore_same_strain_cooccur_clusters.rds")
Snodgrassella_noncore_same_strain_cooccur_clusters <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Snodgrassella_noncore_same_strain_cooccur_clusters.rds")

names(Gilliamella_noncore_same_strain_cooccur_clusters) <- paste("Gilliamella_cluster_", as.character(1:length(Gilliamella_noncore_same_strain_cooccur_clusters)), sep = "")
names(Snodgrassella_noncore_same_strain_cooccur_clusters) <- paste("Snodgrassella_cluster_", as.character(1:length(Snodgrassella_noncore_same_strain_cooccur_clusters)), sep = "")

Gilliamella_and_Snodgrassella_present_noncore$cluster <- rownames(Gilliamella_and_Snodgrassella_present_noncore)

for (Gilliamella_cluster_id in names(Gilliamella_noncore_same_strain_cooccur_clusters)) {
    Gilliamella_and_Snodgrassella_present_noncore[Gilliamella_noncore_same_strain_cooccur_clusters[[Gilliamella_cluster_id]], "cluster"] <- Gilliamella_cluster_id
}

for (Snodgrassella_cluster_id in names(Snodgrassella_noncore_same_strain_cooccur_clusters)) {
  Gilliamella_and_Snodgrassella_present_noncore[Snodgrassella_noncore_same_strain_cooccur_clusters[[Snodgrassella_cluster_id]], "cluster"] <- Snodgrassella_cluster_id
}

Gilliamella_and_Snodgrassella_present_noncore_cluster <- aggregate(. ~ cluster, data = Gilliamella_and_Snodgrassella_present_noncore, FUN = mean)

rownames(Gilliamella_and_Snodgrassella_present_noncore_cluster) <- Gilliamella_and_Snodgrassella_present_noncore_cluster$cluster
Gilliamella_and_Snodgrassella_present_noncore_cluster <- Gilliamella_and_Snodgrassella_present_noncore_cluster[, -which(colnames(Gilliamella_and_Snodgrassella_present_noncore_cluster) == "cluster")]

saveRDS(object = Gilliamella_and_Snodgrassella_present_noncore_cluster,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_and_Snodgrassella_present_noncore_clustered.rds")
