rm(list = ls(all.names = TRUE))

library("reshape2")
library("ggnetwork")
library("igraph")

discover_output <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_cooccur_results_clustered_DBH0.3.rds")

Bifidobacterium_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Bifidobacterium/gene_presence_absence_roary.csv"

Bifidobacterium_panaroo_out <- read.table(Bifidobacterium_path_to_panaroo,
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Bifidobacterium_panaroo_out) <- paste("Bifidobacterium", rownames(Bifidobacterium_panaroo_out), sep = "_")


Firm4_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Firm4/gene_presence_absence_roary.csv.gz"

Firm4_panaroo_out <- read.table(Firm4_path_to_panaroo,
                                header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Firm4_panaroo_out) <- paste("Firm4", rownames(Firm4_panaroo_out), sep = "_")


Firm5_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Firm5/gene_presence_absence_roary.csv"

Firm5_panaroo_out <- read.table(Firm5_path_to_panaroo,
                                header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Firm5_panaroo_out) <- paste("Firm5", rownames(Firm5_panaroo_out), sep = "_")

Gilliamella_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Gilliamella/gene_presence_absence_roary.csv"
Gilliamella_panaroo_out <- read.table(Gilliamella_path_to_panaroo,
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilliamella", rownames(Gilliamella_panaroo_out), sep = "_")

Snodgrassella_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Snodgrassella/gene_presence_absence_roary.csv"
Snodgrassella_panaroo_out <- read.table(Snodgrassella_path_to_panaroo,
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snodgrassella", rownames(Snodgrassella_panaroo_out), sep = "_")

panaroo_out <- rbind(Bifidobacterium_panaroo_out[, c(1:13)],
                     Firm4_panaroo_out[, c(1:13)],
                     Firm5_panaroo_out[, c(1:13)],
                     Gilliamella_panaroo_out[, c(1:13)],
                     Snodgrassella_panaroo_out[, c(1:13)])

num_samples <- 71

discover_output_edges <- discover_output[, c("gene1", "gene2")]
colnames(discover_output_edges) <- c("from", "to")

all_unique_nodes <- c(discover_output_edges$from, discover_output_edges$to) 
all_unique_nodes <- all_unique_nodes[-which(duplicated(all_unique_nodes))]

nodes_info <- data.frame(id = all_unique_nodes,
                         annot = panaroo_out[all_unique_nodes, "Annotation"])

nodes_info$genus <- NA
nodes_info[grep("Bifidobacterium", nodes_info$id), "genus"] <- "Bifidobacterium"
nodes_info[grep("Firm4", nodes_info$id), "genus"] <- "Firm4"
nodes_info[grep("Firm5", nodes_info$id), "genus"] <- "Firm5"
nodes_info[grep("Gilliamella", nodes_info$id), "genus"] <- "Gilliamella"
nodes_info[grep("Snodgrassella", nodes_info$id), "genus"] <- "Snodgrassella"

net <- graph_from_data_frame(d = discover_output_edges,
                             vertices = nodes_info,
                             directed = FALSE)

saveRDS(object = net,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_cooccur_results_clustered_DBH0.3_net.rds")

# Compute node centrality metrics
net_degree <- degree(net)
net_betweenness <- betweenness(net)

head(panaroo_out[names(sort(net_degree, decreasing = TRUE)), "Annotation"], 10)
head(panaroo_out[names(sort(net_betweenness, decreasing = TRUE)), "Annotation"], 10)

#net_closeness <- closeness(net) # SKIPPED BECAUSE IT DOESN'T MAKE SENSE FOR DISCONNECTED GRAPHS

# Skipping eigen centrality too as it results in very similar results to just looking at degrees (at least based on my test)
# net_eigen_centrality <- eigen_centrality(net)$vector

# Get network components

net_components <- components(net)

saveRDS(object = net_components,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_cooccur_results_clustered_DBH0.3_components.rds")

net_flattened <- ggnetwork(net)

ggplot(net_flattened, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges() +
  geom_nodes(aes(color = genus), size = 3) +
  theme_blank()


# Get total number of genes involved
total_genes <- c(discover_output$gene1, discover_output$gene2)
total_genes <- total_genes[-which(duplicated(total_genes))]

### Look at cross-phylotype cooccurences specifically.
discover_output$phylotype1 <- gsub("_.*", "",discover_output$gene1)
discover_output$phylotype2 <- gsub("_.*", "",discover_output$gene2)



discover_output_cross <- discover_output[which(discover_output$phylotype1 != discover_output$phylotype2), ]

total_cross_genes <- c(discover_output_cross$gene1, discover_output_cross$gene2)
total_cross_genes <- total_cross_genes[-which(duplicated(total_cross_genes))]

cross_gene_annot <- panaroo_out[total_cross_genes, "Annotation"]
cross_gene_annot_defined <- cross_gene_annot[-which(cross_gene_annot == "hypothetical protein")]

discover_output_cross$annot1 <- panaroo_out[discover_output_cross$gene1, "Annotation"]
discover_output_cross$annot2 <- panaroo_out[discover_output_cross$gene2, "Annotation"]

discover_output_cross_defined <- discover_output_cross[-which(discover_output_cross$annot1 == "hypothetical protein"), ]
discover_output_cross_defined <- discover_output_cross_defined[-which(discover_output_cross_defined$annot2 == "hypothetical protein"), ]

write.table(x = discover_output_cross_defined,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_cooccur_results_clustered_DBH0.3_cross_phylotype_defined.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)


# Numbers of cross-phylotype cooccurences of each type.
cross_phylotype_instances <- c()

for (i in 1:nrow(discover_output_cross)) {
  phylotypes <- sort(c(discover_output_cross$phylotype1[i], discover_output_cross$phylotype2[i]))
  cross_phylotype_instances <- c(cross_phylotype_instances, paste(phylotypes, collapse = "_"))
}

tested_genes <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_tested_genes.rds")
tested_genes <- names(tested_genes[which(tested_genes)])

network_genes_only <- total_genes

compute_unnorm_expected_prop <- function(tested_genes, phylotype1, phylotype2) {

    phylotype1_prop <- length(grep(phylotype1, tested_genes)) / length(tested_genes)
    phylotype2_prop <- length(grep(phylotype2, tested_genes)) / length(tested_genes)
  
    return(phylotype1_prop * phylotype2_prop)
}

crossphylo_summary <- data.frame(matrix(NA, nrow = 10, ncol = 4))
colnames(crossphylo_summary) <- c("phylo_group", "observed", "unnorm_expected_prop", "unnorm_expected_prop_sig_genes")

phylotypes <- c("Bifidobacterium", "Firm4", "Firm5", "Gilliamella", "Snodgrassella")

row_i <- 1

for (p1 in phylotypes) {
 
  other_phylotypes <- phylotypes[-c(1:which(phylotypes == p1))]
  
  if (length(other_phylotypes) == 0) { break }
  
  for (p2 in other_phylotypes) {
   
    crossphylo_summary[row_i, "phylo_group"] <- paste(c(p1, p2), collapse = " - ")
    
    crossphylo_summary[row_i, "observed"] <- length(which(cross_phylotype_instances == paste(c(p1, p2), collapse = "_")))

    crossphylo_summary[row_i, "unnorm_expected_prop_all_genes"] <- compute_unnorm_expected_prop(tested_genes = tested_genes,
                                                                                                phylotype1 = p1,
                                                                                                phylotype2 = p2)

    crossphylo_summary[row_i, "unnorm_expected_prop_sig_genes"] <- compute_unnorm_expected_prop(tested_genes = network_genes_only,
                                                                                                phylotype1 = p1,
                                                                                                phylotype2 = p2)
    
    row_i <- row_i + 1          
  }
}

crossphylo_summary$expected_all_prop <- crossphylo_summary$unnorm_expected_prop_all_genes / sum(crossphylo_summary$unnorm_expected_prop_all_genes)
crossphylo_summary$expected_all_genes <- crossphylo_summary$expected_all_prop * nrow(discover_output_cross)

crossphylo_summary$expected_sig_prop <- crossphylo_summary$unnorm_expected_prop_sig_genes / sum(crossphylo_summary$unnorm_expected_prop_sig_genes)
crossphylo_summary$expected_sig_genes <- crossphylo_summary$expected_sig_prop * nrow(discover_output_cross)

crossphylo_summary_melt <- melt(crossphylo_summary[, c("phylo_group", "observed", "expected_all_genes", "expected_sig_genes")], id = "phylo_group", variable.name = "count_type")

ggplot(data = crossphylo_summary_melt, aes(y = phylo_group, x = value, colour = count_type)) +
  geom_point(size = 4) +
  xlab("Number of associations") +
  scale_y_discrete(limits = rev(levels(as.factor(crossphylo_summary_melt$phylo_group)))) +
  scale_colour_manual(values = c("green", "black", "brown"))
