rm(list = ls(all.names = TRUE))

# Scan for cases where alleles have high inter-sample divergence, but where the # of alleles per sample is low.

library("phyloseq")
library("parallel")

StrainFinder_core_gene_info <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/StrainFinder_core_genes_and_samples_inputs.rds")

all_genes_abun <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/all_genes_OTU.rds")

gene_presence_noncore <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.matrix", header = TRUE, sep = "\t", row.names = 1)

species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species_present.txt", stringsAsFactors = FALSE, header = FALSE)$V1

compute_intersample_allelic_dists <- function(g, sp) {
  
  gene_matrix <- all_genes_abun[[sp]][[g]]
  
  intersecting_samples <- which(rownames(gene_matrix) %in% StrainFinder_core_gene_info[[sp]]$samples)
  
  if (length(intersecting_samples) < 2 | ncol(gene_matrix) < 2) { return(NULL) }
  
  gene_matrix <- gene_matrix[intersecting_samples, ]
  
  gene_matrix_binary <- gene_matrix
  gene_matrix_binary[gene_matrix_binary > 0] <- 1
  
  gene_matrix_binary_otu <- otu_table(gene_matrix_binary, taxa_are_rows = FALSE)
  gene_matrix_otu <- otu_table(gene_matrix, taxa_are_rows = FALSE)
  
  jaccard_dist <- as.matrix(distance(gene_matrix_binary_otu, method = "jaccard"))
  simpson_dist <- as.matrix(distance(gene_matrix_binary_otu, method = "sim"))
  jsd_dist <- as.matrix(distance(gene_matrix_otu, method = "jsd"))
  
  raw_dist = list(jaccard = jaccard_dist,
                  simpson = simpson_dist,
                  jsd = jsd_dist)
  
  return(raw_dist)
  
}

all_raw_dist <- list()

for (sp in species) {
  
  print(sp)
  
  sp_noncore_genes <- names(all_genes_abun[[sp]])[which(names(all_genes_abun[[sp]]) %in% rownames(gene_presence_noncore))]
  
  all_raw_dist[[sp]] <- mclapply(sp_noncore_genes, compute_intersample_allelic_dists, sp = sp, mc.cores = 50)
  
  names(all_raw_dist[[sp]]) <- sp_noncore_genes
  
}

raw_sp_gene_info <- list()

for (sp in species) {

  sp_noncore_genes <- names(all_genes_abun[[sp]])[which(names(all_genes_abun[[sp]]) %in% rownames(gene_presence_noncore))]
  
  sp_tmp <- data.frame(matrix(NA, nrow = length(sp_noncore_genes), ncol = 9))
  
  rownames(sp_tmp) <- sp_noncore_genes
  colnames(sp_tmp) <- c("Species", "Gene", "Total_alleles", "Total_samples", "Mean_num_alleles", "Max_num_alleles", "Jaccard_mean", "Simpson_mean", "JSD_mean")
  
  sp_tmp$Species <- sp
  sp_tmp$Gene <- sp_noncore_genes
  
  sp_tmp$Total_alleles <- sapply(sp_noncore_genes, function(x) { ncol(all_genes_abun[[sp]][[x]]) })
  sp_tmp$Total_samples <- sapply(sp_noncore_genes, function(x) { nrow(all_genes_abun[[sp]][[x]]) })
  sp_tmp$Mean_num_alleles <- sapply(sp_noncore_genes, function(x) { mean(rowSums(all_genes_abun[[sp]][[x]] > 0)) })
  sp_tmp$Max_num_alleles <- sapply(sp_noncore_genes, function(x) { max(rowSums(all_genes_abun[[sp]][[x]] > 0)) })
  
  sp_tmp$Jaccard_mean <- sapply(sp_noncore_genes, function(x) { if (is.null(all_raw_dist[[sp]][[x]])) { return(NA) } else { mean(as.dist(all_raw_dist[[sp]][[x]]$jaccard)) } })
  sp_tmp$Simpson_mean <- sapply(sp_noncore_genes, function(x) { if (is.null(all_raw_dist[[sp]][[x]])) { return(NA) } else { mean(as.dist(all_raw_dist[[sp]][[x]]$simpson)) } })
  sp_tmp$JSD_mean <- sapply(sp_noncore_genes, function(x) { if (is.null(all_raw_dist[[sp]][[x]])) { return(NA) } else { mean(as.dist(all_raw_dist[[sp]][[x]]$jsd)) } })
  
  raw_sp_gene_info[[sp]] <- sp_tmp
}

all_gene_info <- do.call(rbind, raw_sp_gene_info)

all_gene_info_10samples_10alleles <- all_gene_info[which(all_gene_info$Total_alleles >= 10 & all_gene_info$Total_alleles >= 10), ]

all_gene_info_10samples_10alleles_Bifidobacterium_asteroides <- all_gene_info_10samples_10alleles[which(all_gene_info_10samples_10alleles$Species == "Bifidobacterium_asteroides"), ]

COG_category_to_gene <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/COG_category_to_ortholog.rds")

all_gene_info_10samples_10alleles$G_category <- "No"
all_gene_info_10samples_10alleles$G_category[which(all_gene_info_10samples_10alleles$Gene %in% COG_category_to_gene$G)] <- "Yes"

ggplot(data = all_gene_info_10samples_10alleles, aes(x = Mean_num_alleles, y = JSD_mean)) +
  geom_point(alpha = 0.1) +
  theme_bw() +
  xlab("Mean number of alleles per sample") +
  ylab("Mean Jensen-Shannon Divergence") +
  coord_cartesian(xlim = c(0, 17), ylim = c(0, 0.6))


ggplot(data = all_gene_info_10samples_10alleles, aes(x = Max_num_alleles, y = JSD_mean)) +
  geom_point(alpha = 0.1) +
  theme_bw() +
  xlab("Max number of alleles per sample") +
  ylab("Mean Jensen-Shannon Divergence") +
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 0.6))

combined_panaroo <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/combined_panaroo_annot.rds")

combined_panaroo["Snodgrassella_alvi_spuE_1", ]

ggplot(data = all_gene_info_10samples_10alleles_Bifidobacterium_asteroides, aes(x = Mean_num_alleles, y = JSD_mean)) +
  geom_point() +
  facet_wrap(Species ~ .)


tmp <- all_gene_info_10samples_10alleles_Bifidobacterium_asteroides[which(all_gene_info_10samples_10alleles_Bifidobacterium_asteroides$JSD_mean >= 0.45 & all_gene_info_10samples_10alleles_Bifidobacterium_asteroides$Mean_num_alleles < 5), ]

ortholog_COG_categories <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/all_microbiota_COG.category_by_ortholog.tsv",
                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
ortholog_COG_categories[tmp$Gene, ]


combined_panaroo[tmp$Gene, "Annotation"]