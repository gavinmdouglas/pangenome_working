rm(list = ls(all.names = TRUE))

library("reshape2")
library("ggnetwork")
library("igraph")

discover_output_cooccur <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_cooccur_results_clustered_DBH0.35.rds")
discover_output_mutex <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_mutex_results_clustered_DBH0.35.rds")

species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species.txt", stringsAsFactors = FALSE)$V1

num_samples <- 74

panaroo_out <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/combined_panaroo_annot.rds")


discover_output_cooccur_edges <- discover_output_cooccur[, c("gene1", "gene2")]
colnames(discover_output_cooccur_edges) <- c("from", "to")

all_unique_nodes <- c(discover_output_cooccur_edges$from, discover_output_cooccur_edges$to) 
all_unique_nodes <- all_unique_nodes[-which(duplicated(all_unique_nodes))]

nodes_info_cooccur <- data.frame(id = all_unique_nodes,
                         annot = panaroo_out[all_unique_nodes, "Annotation"])

nodes_info_cooccur$species <- NA

for (sp in species) {
 sp_matches <-  grep(sp, nodes_info_cooccur$id)
 
 if (length(sp_matches) > 0) {
   nodes_info_cooccur[sp_matches, "species"] <- sp
 }
}

net_cooccur <- graph_from_data_frame(d = discover_output_cooccur_edges,
                             vertices = nodes_info_cooccur,
                             directed = FALSE)

saveRDS(object = net_cooccur,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_cooccur_results_clustered_DBH0.35_net.rds")

# Compute node centrality metrics
net_cooccur_degree <- degree(net_cooccur)
net_cooccur_betweenness <- betweenness(net_cooccur)

head(panaroo_out[names(sort(net_cooccur_degree, decreasing = TRUE)), "Annotation"], 10)
head(panaroo_out[names(sort(net_cooccur_betweenness, decreasing = TRUE)), "Annotation"], 10)

#net_cooccur_closeness <- closeness(net_cooccur) # SKIPPED BECAUSE IT DOESN'T MAKE SENSE FOR DISCONNECTED GRAPHS

# Skipping eigen centrality too as it results in very similar results to just looking at degrees (at least based on my test)
# net_cooccur_eigen_centrality <- eigen_centrality(net_cooccur)$vector

# Get net_cooccurwork components

net_cooccur_components <- components(net_cooccur)

saveRDS(object = net_cooccur_components,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_cooccur_results_clustered_DBH0.35_net_components.rds")

net_cooccur_flattened <- ggnetwork(net_cooccur)

ggplot(net_cooccur_flattened, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges() +
  geom_nodes(aes(color = species), size = 3) +
  theme_blank()


# Get total number of genes involved
total_genes_cooccur <- c(discover_output_cooccur$gene1, discover_output_cooccur$gene2)
total_genes_cooccur <- total_genes_cooccur[-which(duplicated(total_genes_cooccur))]

### Look at cross-species cooccurences specifically.
discover_output_cooccur$species1 <- NA
discover_output_cooccur$species2 <- NA

for (sp in species) {
  
  if (length(grep(sp, discover_output_cooccur$gene1)) > 0) {
    discover_output_cooccur[grep(sp, discover_output_cooccur$gene1), "species1"] <- sp
  }
  
  if (length(grep(sp, discover_output_cooccur$gene2)) > 0) {
    discover_output_cooccur[grep(sp, discover_output_cooccur$gene2), "species2"] <- sp
  }
}



discover_output_cooccur_cross <- discover_output_cooccur[which(discover_output_cooccur$species1 != discover_output_cooccur$species2), ]

total_cross_genes <- c(discover_output_cooccur_cross$gene1, discover_output_cooccur_cross$gene2)
total_cross_genes <- total_cross_genes[-which(duplicated(total_cross_genes))]

cross_gene_annot <- panaroo_out[total_cross_genes, "Annotation"]

discover_output_cooccur_cross$annot1 <- panaroo_out[discover_output_cooccur_cross$gene1, "Annotation"]
discover_output_cooccur_cross$annot2 <- panaroo_out[discover_output_cooccur_cross$gene2, "Annotation"]

write.table(x = discover_output_cooccur_cross,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_cooccur_results_clustered_DBH0.35_cross_species.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)


# Numbers of cross-species cooccurences of each type.
cross_species_instances <- c()

for (i in 1:nrow(discover_output_cooccur_cross)) {
  species_rep <- sort(c(discover_output_cooccur_cross$species1[i], discover_output_cooccur_cross$species2[i]))
  cross_species_instances <- c(cross_species_instances, paste(species_rep, collapse = "_"))
}

tested_genes <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_tested_genes.rds")
tested_genes <- names(tested_genes[which(tested_genes)])

network_genes_only <- total_genes_cooccur

compute_unnorm_expected_prop <- function(tested_genes, species1, species2) {

    species1_prop <- length(grep(species1, tested_genes)) / length(tested_genes)
    species2_prop <- length(grep(species2, tested_genes)) / length(tested_genes)
  
    return(species1_prop * species2_prop)
}

crossphylo_summary <- data.frame(matrix(NA, nrow = 0, ncol = 4))
colnames(crossphylo_summary) <- c("species_group", "observed", "unnorm_expected_prop", "unnorm_expected_prop_sig_genes")

row_i <- 1

species_present <- unique(c(discover_output_cooccur_cross$species1, discover_output_cooccur_cross$species2))

for (p1 in species_present) {
 
  other_species <- species_present[-c(1:which(species_present == p1))]
  
  if (length(other_species) == 0) { break }
  
  for (p2 in other_species) {
   
    crossphylo_summary[row_i, "species_group"] <- paste(c(p1, p2), collapse = " - ")
    
    crossphylo_summary[row_i, "observed"] <- length(which(cross_species_instances == paste(c(p1, p2), collapse = "_")))

    crossphylo_summary[row_i, "unnorm_expected_prop_all_genes"] <- compute_unnorm_expected_prop(tested_genes = tested_genes,
                                                                                                species1 = p1,
                                                                                                species2 = p2)

    crossphylo_summary[row_i, "unnorm_expected_prop_sig_genes"] <- compute_unnorm_expected_prop(tested_genes = network_genes_only,
                                                                                                species1 = p1,
                                                                                                species2 = p2)
    
    row_i <- row_i + 1          
  }
}

crossphylo_summary$expected_all_prop <- crossphylo_summary$unnorm_expected_prop_all_genes / sum(crossphylo_summary$unnorm_expected_prop_all_genes)
crossphylo_summary$expected_all_genes <- crossphylo_summary$expected_all_prop * nrow(discover_output_cooccur_cross)

crossphylo_summary$expected_sig_prop <- crossphylo_summary$unnorm_expected_prop_sig_genes / sum(crossphylo_summary$unnorm_expected_prop_sig_genes)
crossphylo_summary$expected_sig_genes <- crossphylo_summary$expected_sig_prop * nrow(discover_output_cooccur_cross)

crossphylo_summary_melt <- melt(crossphylo_summary[, c("species_group", "observed", "expected_all_genes", "expected_sig_genes")], id = "species_group", variable.name = "count_type")

ggplot(data = crossphylo_summary_melt, aes(y = species_group, x = value, colour = count_type)) +
  geom_point(size = 4) +
  xlab("Number of associations") +
  scale_y_discrete(limits = rev(levels(as.factor(crossphylo_summary_melt$species_group)))) +
  scale_colour_manual(values = c("green", "black", "brown"))






### Get key files for mutually exclusive features too

discover_output_mutex_edges <- discover_output_mutex[, c("gene1", "gene2")]
colnames(discover_output_mutex_edges) <- c("from", "to")

all_unique_nodes <- c(discover_output_mutex_edges$from, discover_output_mutex_edges$to) 
all_unique_nodes <- all_unique_nodes[-which(duplicated(all_unique_nodes))]

nodes_info_mutex <- data.frame(id = all_unique_nodes,
                                 annot = panaroo_out[all_unique_nodes, "Annotation"])

nodes_info_mutex$species <- NA

for (sp in species) {
  sp_matches <-  grep(sp, nodes_info_mutex$id)
  
  if (length(sp_matches) > 0) {
    nodes_info_mutex[sp_matches, "species"] <- sp
  }
}

net_mutex <- graph_from_data_frame(d = discover_output_mutex_edges,
                                     vertices = nodes_info_mutex,
                                     directed = FALSE)

saveRDS(object = net_mutex,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_mutex_results_clustered_DBH0.35_net.rds")

# Compute node centrality metrics
net_mutex_degree <- degree(net_mutex)
net_mutex_betweenness <- betweenness(net_mutex)

head(panaroo_out[names(sort(net_mutex_degree, decreasing = TRUE)), "Annotation"], 10)
head(panaroo_out[names(sort(net_mutex_betweenness, decreasing = TRUE)), "Annotation"], 10)

#net_mutex_closeness <- closeness(net_mutex) # SKIPPED BECAUSE IT DOESN'T MAKE SENSE FOR DISCONNECTED GRAPHS

# Skipping eigen centrality too as it results in very similar results to just looking at degrees (at least based on my test)
# net_mutex_eigen_centrality <- eigen_centrality(net_mutex)$vector

# Get net_mutexwork components

net_mutex_components <- components(net_mutex)

saveRDS(object = net_mutex_components,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_mutex_results_clustered_DBH0.35_net_components.rds")

net_mutex_flattened <- ggnetwork(net_mutex)

ggplot(net_mutex_flattened, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges() +
  geom_nodes(aes(color = species), size = 3) +
  theme_blank()


total_genes_mutex <- c(discover_output_mutex$gene1, discover_output_mutex$gene2)
total_genes_mutex <- total_genes_mutex[-which(duplicated(total_genes_mutex))]

### Look at cross-species mutexences specifically.
discover_output_mutex$species1 <- NA
discover_output_mutex$species2 <- NA

for (sp in species) {
  
  if (length(grep(sp, discover_output_mutex$gene1)) > 0) {
    discover_output_mutex[grep(sp, discover_output_mutex$gene1), "species1"] <- sp
  }
  
  if (length(grep(sp, discover_output_mutex$gene2)) > 0) {
    discover_output_mutex[grep(sp, discover_output_mutex$gene2), "species2"] <- sp
  }
}



discover_output_mutex_cross <- discover_output_mutex[which(discover_output_mutex$species1 != discover_output_mutex$species2), ]

total_cross_genes <- c(discover_output_mutex_cross$gene1, discover_output_mutex_cross$gene2)
total_cross_genes <- total_cross_genes[-which(duplicated(total_cross_genes))]

cross_gene_annot <- panaroo_out[total_cross_genes, "Annotation"]

discover_output_mutex_cross$annot1 <- panaroo_out[discover_output_mutex_cross$gene1, "Annotation"]
discover_output_mutex_cross$annot2 <- panaroo_out[discover_output_mutex_cross$gene2, "Annotation"]

write.table(x = discover_output_mutex_cross,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_mutex_results_clustered_DBH0.35_cross_species.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)


# Numbers of cross-species mutexences of each type.
cross_species_instances <- c()

for (i in 1:nrow(discover_output_mutex_cross)) {
  species_rep <- sort(c(discover_output_mutex_cross$species1[i], discover_output_mutex_cross$species2[i]))
  cross_species_instances <- c(cross_species_instances, paste(species_rep, collapse = "_"))
}

tested_genes <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_tested_genes.rds")

network_genes_only <- total_genes_mutex

compute_unnorm_expected_prop <- function(tested_genes, species1, species2) {
  
  species1_prop <- length(grep(species1, tested_genes)) / length(tested_genes)
  species2_prop <- length(grep(species2, tested_genes)) / length(tested_genes)
  
  return(species1_prop * species2_prop)
}

crossphylo_summary <- data.frame(matrix(NA, nrow = 0, ncol = 4))
colnames(crossphylo_summary) <- c("species_group", "observed", "unnorm_expected_prop", "unnorm_expected_prop_sig_genes")

row_i <- 1

species_present <- unique(c(discover_output_mutex_cross$species1, discover_output_mutex_cross$species2))

for (p1 in species_present) {
  
  other_species <- species_present[-c(1:which(species_present == p1))]
  
  if (length(other_species) == 0) { break }
  
  for (p2 in other_species) {
    
    crossphylo_summary[row_i, "species_group"] <- paste(c(p1, p2), collapse = " - ")
    
    crossphylo_summary[row_i, "observed"] <- length(which(cross_species_instances == paste(c(p1, p2), collapse = "_")))
    
    crossphylo_summary[row_i, "unnorm_expected_prop_all_genes"] <- compute_unnorm_expected_prop(tested_genes = tested_genes,
                                                                                                species1 = p1,
                                                                                                species2 = p2)
    
    crossphylo_summary[row_i, "unnorm_expected_prop_sig_genes"] <- compute_unnorm_expected_prop(tested_genes = network_genes_only,
                                                                                                species1 = p1,
                                                                                                species2 = p2)
    
    row_i <- row_i + 1          
  }
}

crossphylo_summary$expected_all_prop <- crossphylo_summary$unnorm_expected_prop_all_genes / sum(crossphylo_summary$unnorm_expected_prop_all_genes)
crossphylo_summary$expected_all_genes <- crossphylo_summary$expected_all_prop * nrow(discover_output_mutex_cross)

crossphylo_summary$expected_sig_prop <- crossphylo_summary$unnorm_expected_prop_sig_genes / sum(crossphylo_summary$unnorm_expected_prop_sig_genes)
crossphylo_summary$expected_sig_genes <- crossphylo_summary$expected_sig_prop * nrow(discover_output_mutex_cross)

crossphylo_summary_melt <- melt(crossphylo_summary[, c("species_group", "observed", "expected_all_genes", "expected_sig_genes")], id = "species_group", variable.name = "count_type")

ggplot(data = crossphylo_summary_melt, aes(y = species_group, x = value, colour = count_type)) +
  geom_point(size = 4) +
  xlab("Number of associations") +
  scale_y_discrete(limits = rev(levels(as.factor(crossphylo_summary_melt$species_group)))) +
  scale_colour_manual(values = c("green", "black", "brown"))

