rm(list = ls(all.names = TRUE))

library("ggnetwork")
library("igraph")
library("RCy3")

discover_output <- readRDS("/Users/Gavin/R_projects/honey_bee_pangenome/data/Ellegaard_pangenome_coincidence/Ellegaard_2019_2020_DISCOVER_cooccur_results_clustered_DBH0.3.rds")

Gilliamella_panaroo_out <- read.table("honey_bee_pangenome/data/binning_output/panaroo_output_dastools_and_ref/Gilliamella_gene_presence_absence_roary-formatted.csv",
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilli", rownames(Gilliamella_panaroo_out), sep = "_")

Snodgrassella_panaroo_out <- read.table("honey_bee_pangenome/data/binning_output/panaroo_output_dastools_and_ref/Snodgrassella_gene_presence_absence_roary-formatted.csv",
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snod", rownames(Snodgrassella_panaroo_out), sep = "_")

panaroo_out <- rbind(Gilliamella_panaroo_out[, c(1:13)], Snodgrassella_panaroo_out[, c(1:13)])

num_samples <- 74

discover_output_edges <- discover_output[, c("gene1", "gene2")]
colnames(discover_output_edges) <- c("from", "to")

all_unique_nodes <- c(discover_output_edges$from, discover_output_edges$to) 
all_unique_nodes <- all_unique_nodes[-which(duplicated(all_unique_nodes))]

nodes_info <- data.frame(id = all_unique_nodes,
                         annot = panaroo_out[all_unique_nodes, "Annotation"])

nodes_info$genus <- NA
nodes_info[grep("Gilli", nodes_info$id), "genus"] <- "Gilliamella"
nodes_info[grep("Snod", nodes_info$id), "genus"] <- "Snodgrassella"

net <- graph_from_data_frame(d = discover_output_edges,
                             vertices = nodes_info,
                             directed = FALSE)

nodes_genus_col <- rep(NA, nrow(nodes_info))
nodes_genus_col[which(nodes_info$genus == "Gilliamella")] <- "blue"
nodes_genus_col[which(nodes_info$genus == "Snodgrassella")] <- "green"

# Quick way to plot with igraph
#plot(net, vertex.color = nodes_genus_col)

# Compute node centrality metrics
net_degree <- degree(net)
net_betweenness <- betweenness(net)

head(panaroo_out[names(sort(net_degree, decreasing = TRUE)), "Annotation"], 20)
head(panaroo_out[names(sort(net_betweenness, decreasing = TRUE)), "Annotation"], 20)

#net_closeness <- closeness(net) # SKIPPED BECAUSE IT DOESN'T MAKE SENSE FOR DISCONNECTED GRAPHS

# Skipping eigen centrality too as it results in very similar results to just looking at degrees (at least based on my test)
# net_eigen_centrality <- eigen_centrality(net)$vector

# Get network components

components(net)

net_flattened <- ggnetwork(net)

ggplot(net_flattened, aes(x = x, y = y, xend = xend, yend = yend)) +
       geom_edges() +
       geom_nodes(aes(color = genus), size = 1) +
       theme_blank()

# To open graph in Cytoscape, which could be useful for figuring out the functional similarities of genes in the same cluster
createNetworkFromIgraph(net, "myGraph")
