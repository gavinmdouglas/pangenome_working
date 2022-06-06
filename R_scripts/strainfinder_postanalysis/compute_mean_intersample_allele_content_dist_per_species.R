rm(list = ls(all.names = TRUE))

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



num_comparisons_between_samples <- list()

summed_jaccard_dist_between_samples <- list()
summed_simpson_dist_between_samples <- list()
summed_jsd_dist_between_samples <- list()

mean_jaccard_dist <- list()
mean_simpson_dist <- list()
mean_jsd_dist <- list()

for (sp in species) {
 
  print(sp)
  
  empty_df <- data.frame(matrix(0,
                                nrow = length(StrainFinder_core_gene_info[[sp]]$samples),
                                ncol = length(StrainFinder_core_gene_info[[sp]]$samples)))
  
  rownames(empty_df) <- StrainFinder_core_gene_info[[sp]]$samples
  colnames(empty_df) <- StrainFinder_core_gene_info[[sp]]$samples
 
  num_comparisons_between_samples[[sp]] <- empty_df

  summed_jaccard_dist_between_samples[[sp]] <- empty_df
  summed_simpson_dist_between_samples[[sp]] <- empty_df
  summed_jsd_dist_between_samples[[sp]] <- empty_df

  sp_noncore_genes <- names(all_genes_abun[[sp]])[which(names(all_genes_abun[[sp]]) %in% rownames(gene_presence_noncore))]
  
  for (g in sp_noncore_genes) {
  
    gene_matrix <- all_genes_abun[[sp]][[g]]
    
    intersecting_samples <- rownames(gene_matrix)[which(rownames(gene_matrix) %in% StrainFinder_core_gene_info[[sp]]$samples)]

    if (length(intersecting_samples) < 2 | ncol(gene_matrix) < 2) { next }
    
    
    jaccard_dist <- all_raw_dist[[sp]][[g]]$jaccard
    simpson_dist <- all_raw_dist[[sp]][[g]]$simpson
    jsd_dist <- all_raw_dist[[sp]][[g]]$jsd
    
    num_comparison_dummy <- jaccard_dist
    num_comparison_dummy[num_comparison_dummy > 0] <- 1
    
    
    
    summed_jaccard_dist_between_samples[[sp]][intersecting_samples, intersecting_samples] <- summed_jaccard_dist_between_samples[[sp]][intersecting_samples, intersecting_samples] + jaccard_dist
    summed_simpson_dist_between_samples[[sp]][intersecting_samples, intersecting_samples] <- summed_simpson_dist_between_samples[[sp]][intersecting_samples, intersecting_samples] + simpson_dist
    summed_jsd_dist_between_samples[[sp]][intersecting_samples, intersecting_samples] <- summed_jsd_dist_between_samples[[sp]][intersecting_samples, intersecting_samples] + jsd_dist
    num_comparisons_between_samples[[sp]][intersecting_samples, intersecting_samples] <- num_comparisons_between_samples[[sp]][intersecting_samples, intersecting_samples] + num_comparison_dummy
  }
  
  mean_jaccard_dist[[sp]] <- as.dist(summed_jaccard_dist_between_samples[[sp]] / num_comparisons_between_samples[[sp]])
  mean_simpson_dist[[sp]] <- as.dist(summed_simpson_dist_between_samples[[sp]] / num_comparisons_between_samples[[sp]])
  mean_jsd_dist[[sp]] <- as.dist(summed_jsd_dist_between_samples[[sp]] / num_comparisons_between_samples[[sp]])
  
}


mean_allelic_dists <- list(jaccard = mean_jaccard_dist,
                           simpson = mean_simpson_dist,
                           jsd = mean_jsd_dist)

saveRDS(object = mean_allelic_dists,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/allelic_mean_intersample_dists.rds")
