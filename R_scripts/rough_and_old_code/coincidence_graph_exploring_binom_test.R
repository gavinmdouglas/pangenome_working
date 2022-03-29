rm(list = ls(all.names = TRUE))

library("ggnetwork")
library("igraph")
library("RCy3")

binom_test_output <- readRDS("/Users/Gavin/R_projects/honey_bee_pangenome/data/Ellegaard_pangenome_coincidence/Gilliamella_and_Snodgrassella_panaroo_coincidence_binomial_BH_0.4_df.rds")

Gilliamella_panaroo_out <- read.table("honey_bee_pangenome/data/binning_output/panaroo_output_dastools_and_ref/Gilliamella_gene_presence_absence_roary-formatted.csv",
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilli", rownames(Gilliamella_panaroo_out), sep = "_")

Snodgrassella_panaroo_out <- read.table("honey_bee_pangenome/data/binning_output/panaroo_output_dastools_and_ref/Snodgrassella_gene_presence_absence_roary-formatted.csv",
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snod", rownames(Snodgrassella_panaroo_out), sep = "_")

panaroo_out <- rbind(Gilliamella_panaroo_out[, c(1:13)], Snodgrassella_panaroo_out[, c(1:13)])


binom_test_output$weight <- abs(log((binom_test_output$obs_prop * 54 + 1) / (binom_test_output$exp_prop * 54 + 1)))
# binom_test_output$type <- NA
# binom_test_output[which(binom_test_output$exp_prop > binom_test_output$obs_prop), "type"] <- "negative"
# binom_test_output[which(binom_test_output$exp_prop < binom_test_output$obs_prop), "type"] <- "positive"

binom_test_output_sig <- binom_test_output[which(binom_test_output$BH < 0.4), ]


# Get numbers of intra and inter genus coincidences
Gilli_gene1 <- grep("Gilli", binom_test_output$gene1)
Gilli_gene2 <- grep("Gilli", binom_test_output$gene2)
Snod_gene1 <- grep("Snod", binom_test_output$gene1)
Snod_gene2 <- grep("Snod", binom_test_output$gene2)

length(which(Gilli_gene1 %in% Snod_gene2)) + length(which(Gilli_gene2 %in% Snod_gene1))

binom_test_edges <- binom_test_output_sig[, c("gene1", "gene2", "weight")]
colnames(binom_test_edges) <- c("from", "to", "weight")

all_unique_nodes <- c(binom_test_edges$from, binom_test_edges$to) 
all_unique_nodes <- all_unique_nodes[-which(duplicated(all_unique_nodes))]

nodes_info <- data.frame(id = all_unique_nodes,
                         annot = panaroo_out[all_unique_nodes, "Annotation"])

nodes_info$genus <- NA
nodes_info[grep("Gilli", nodes_info$id), "genus"] <- "Gilliamella"
nodes_info[grep("Snod", nodes_info$id), "genus"] <- "Snodgrassella"

net <- graph_from_data_frame(d = binom_test_edges,
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

#net_closeness <- closeness(net) # SKIPPED BECAUSE IT DOESN'T MAKE SENSE FOR DISCONNECTED GRAPHS

# Skipping eigen centrality too as it results in very similar results to just looking at degrees (at least based on my test)
# net_eigen_centrality <- eigen_centrality(net)$vector

net_flattened <- ggnetwork(net)

net_flattened$weight_rounded <- as.character(round(net_flattened$weight, digits = 1))
unique_weight_rounded <- table(net_flattened$weight_rounded)
weight_col_range <- colorRampPalette(c("grey30", "black"))(length(unique_weight_rounded))

net_flattened$weight_rounded_col <- NA

for (i in 1:length(unique_weight_rounded)) {
  rounded_weight <- names(unique_weight_rounded)[i]
  net_flattened[which(net_flattened$weight_rounded == rounded_weight), "weight_rounded_col"] <- weight_col_range[i]
}
  
ggplot(net_flattened, aes(x = x, y = y, xend = xend, yend = yend)) +
       geom_edges() +
       geom_nodes(aes(color = genus), size = 5) +
       theme_blank()

# To open graph in Cytoscape, which could be useful for figuring out the functional similarities of genes in the same cluster
createNetworkFromIgraph(net, "myGraph")
