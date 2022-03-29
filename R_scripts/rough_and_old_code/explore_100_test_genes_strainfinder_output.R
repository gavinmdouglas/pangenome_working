### Compare number of haplotypes within each sample and how the haplotypes are distributed across samples.

rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(ggplot2)

shannon <- function(x) { 
  x <- x / sum(x)
  return(-1 * sum(x * log(x)))
}

simpson <- function(x) {
  x <- x / sum(x)
  return(1 - sum(x ** 2))
}

invsimpson <- function(x) {
  x <- x / sum(x)
  return(1 / sum(x ** 2))
}

all_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2019_and_2020_SRR_ids.txt",
                          header = FALSE, stringsAsFactors = FALSE)$V1

core_gene_output <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_genes_output/strainfinder_pandora_struct_abun.rds")

panaroo_summary <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/combined_panaroo_and_core_genes.rds")
panaroo_summary$panaroo_all_phylotypes$percent_isolates <- NA

panaroo_num_genomes <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/panaroo_num_ref_genomes.rds")

for (phylotype in names(panaroo_num_genomes)) {
  phylotype_genes <- grep(phylotype, rownames(panaroo_summary$panaroo_all_phylotypes), value = TRUE)
  panaroo_summary$panaroo_all_phylotypes[phylotype_genes, "percent_isolates"] <- (panaroo_summary$panaroo_all_phylotypes[phylotype_genes, "No..isolates"] / panaroo_num_genomes[[phylotype]]) * 100
}

num_strains_per_sample <- list()
num_strains_per_sample[["Bifidobacterium"]] <- rowSums(core_gene_output$Bifidobacterium > 0)
num_strains_per_sample[["Firm4"]] <- rowSums(core_gene_output$Firm4 > 0)
num_strains_per_sample[["Firm5"]] <- rowSums(core_gene_output$Firm5 > 0)
num_strains_per_sample[["Gilliamella"]] <- rowSums(core_gene_output$Gilliamella > 0)
num_strains_per_sample[["Snodgrassella"]] <- rowSums(core_gene_output$Snodgrassella > 0)

strain_shannon_per_sample <- list()
strain_shannon_per_sample[["Bifidobacterium"]] <- apply(core_gene_output$Bifidobacterium, 1, function(x) { x <- x[which(x > 0)]; shannon(x) })
strain_shannon_per_sample[["Firm4"]] <- apply(core_gene_output$Firm4, 1, function(x) { x <- x[which(x > 0)]; shannon(x) })
strain_shannon_per_sample[["Firm5"]] <- apply(core_gene_output$Firm5, 1, function(x) { x <- x[which(x > 0)]; shannon(x) })
strain_shannon_per_sample[["Gilliamella"]] <- apply(core_gene_output$Gilliamella, 1, function(x) { x <- x[which(x > 0)]; shannon(x) })
strain_shannon_per_sample[["Snodgrassella"]] <- apply(core_gene_output$Snodgrassella, 1, function(x) { x <- x[which(x > 0)]; shannon(x) })

strain_simpson_per_sample <- list()
strain_simpson_per_sample[["Bifidobacterium"]] <- apply(core_gene_output$Bifidobacterium, 1, function(x) { x <- x[which(x > 0)]; simpson(x) })
strain_simpson_per_sample[["Firm4"]] <- apply(core_gene_output$Firm4, 1, function(x) { x <- x[which(x > 0)]; simpson(x) })
strain_simpson_per_sample[["Firm5"]] <- apply(core_gene_output$Firm5, 1, function(x) { x <- x[which(x > 0)]; simpson(x) })
strain_simpson_per_sample[["Gilliamella"]] <- apply(core_gene_output$Gilliamella, 1, function(x) { x <- x[which(x > 0)]; simpson(x) })
strain_simpson_per_sample[["Snodgrassella"]] <- apply(core_gene_output$Snodgrassella, 1, function(x) { x <- x[which(x > 0)]; simpson(x) })

strain_invsimpson_per_sample <- list()
strain_invsimpson_per_sample[["Bifidobacterium"]] <- apply(core_gene_output$Bifidobacterium, 1, function(x) { x <- x[which(x > 0)]; invsimpson(x) })
strain_invsimpson_per_sample[["Firm4"]] <- apply(core_gene_output$Firm4, 1, function(x) { x <- x[which(x > 0)]; invsimpson(x) })
strain_invsimpson_per_sample[["Firm5"]] <- apply(core_gene_output$Firm5, 1, function(x) { x <- x[which(x > 0)]; invsimpson(x) })
strain_invsimpson_per_sample[["Gilliamella"]] <- apply(core_gene_output$Gilliamella, 1, function(x) { x <- x[which(x > 0)]; invsimpson(x) })
strain_invsimpson_per_sample[["Snodgrassella"]] <- apply(core_gene_output$Snodgrassella, 1, function(x) { x <- x[which(x > 0)]; invsimpson(x) })


test_genes <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/all_test_genes.txt",
                         stringsAsFactors = FALSE)$V1

data_dir <- "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/"

min_abun <- 0.01

col_categories <- c("phylotype",
                    "sample_presence",
                    "total_num_haplotypes",
                    "min_num_haplotypes", "max_num_haplotypes",
                    "mean_num_haplotypes", "sd_num_haplotypes", "cv_num_haplotypes",
                    "shannon_mean", "shannon_sd", "shannon_cv",
                    "simpson_mean", "simpson_sd", "simpson_cv",
                    "invsimpson_mean", "invsimpson_sd", "invsimpson_cv",
                    "mean_norm_num_haplotypes", "sd_norm_num_haplotypes", "cv_norm_num_haplotypes",
                    "shannon_norm_mean", "shannon_norm_sd", "shannon_norm_cv",
                    "simpson_norm_mean", "simpson_norm_sd", "simpson_norm_cv",
                    "invsimpson_norm_mean", "invsimpson_norm_sd", "invsimpson_norm_cv")


gene_haplotype_metrics <- data.frame(matrix(NA, nrow = length(test_genes), ncol = length(col_categories)))
colnames(gene_haplotype_metrics) <- col_categories
rownames(gene_haplotype_metrics) <- test_genes

for (gene in test_genes) {

  gene_phylotype <- strsplit(gene, split = "_")[[1]][1]
  
  gene_samples_filepath <- paste(data_dir,
                                 "prepped_input/",
                                 gene,
                                 "_samples.txt",
                                 sep = "")
  
  if (file.exists(gene_samples_filepath)) {
  
    gene_samples <- read.table(file = gene_samples_filepath, stringsAsFactors = FALSE)$V1
    
    gene_aic_filepath <- paste(data_dir,
                               "/output/pandora_prepped_w_struct/",
                               gene,
                               ".strain_fit_summary.tsv",
                               sep = "")
    
    gene_aic_output <- read.table(file = gene_aic_filepath, header = TRUE, sep = "\t")
    
    predicted_num_strains <- gene_aic_output[which.min(gene_aic_output$AIC), "ID"]
    
    abun_path <- paste(data_dir,
                       "/output/pandora_prepped_w_struct/otu_tables/otu_table.",
                       gene,
                       ".",
                       as.character(predicted_num_strains),
                       ".txt",
                       sep = "")
    
    abun <- read.table(abun_path, header = FALSE, sep = "\t", skip = 1)
    
    rownames(abun) <- gene_samples
    
    abun_filt <- abun
    
    abun_filt[abun_filt < min_abun] <- 0
    
    abun_filt <- data.frame(sweep(x = abun_filt,
                            MARGIN = 1,
                            STATS = rowSums(abun_filt),
                            FUN = '/')) * 100
    
    if (length(which(colSums(abun_filt) == 0)) > 0) {
      abun_filt <- abun_filt[, -which(colSums(abun_filt) == 0)]
    }
    
    if (length(which(rowSums(abun_filt) == 0)) > 0) {
      abun_filt <- abun_filt[-which(rowSums(abun_filt) == 0), ]
    }
  
  } else {
  
    abun_filt <- data.frame(matrix(1, nrow = length(num_strains_per_sample[[gene_phylotype]]), ncol = 1))
  
  }
  
  num_haplotypes_per_sample <- rowSums(abun_filt > 0)
  gene_shannon_out <- apply(abun_filt, 1, function(x) { x <- x[which(x > 0)]; shannon(x) })
  gene_simpson_out <- apply(abun_filt, 1, function(x) { x <- x[which(x > 0)]; simpson(x) })
  gene_invsimpson_out <- apply(abun_filt, 1, function(x) { x <- x[which(x > 0)]; invsimpson(x) })
  
  samples_w_strains_called <- names(num_haplotypes_per_sample)[which(names(num_haplotypes_per_sample) %in% names(num_strains_per_sample[[gene_phylotype]]))]
  
  norm_num_haplotypes_per_sample <- num_haplotypes_per_sample[samples_w_strains_called] / num_strains_per_sample[[gene_phylotype]][samples_w_strains_called]
  norm_shannon_per_sample <- gene_shannon_out[samples_w_strains_called] / strain_shannon_per_sample[[gene_phylotype]][samples_w_strains_called]
  norm_simpson_per_sample <- gene_simpson_out[samples_w_strains_called] / strain_simpson_per_sample[[gene_phylotype]][samples_w_strains_called]
  norm_invsimpson_per_sample <- gene_invsimpson_out[samples_w_strains_called] / strain_invsimpson_per_sample[[gene_phylotype]][samples_w_strains_called]
  
  gene_haplotype_metrics[gene, "phylotype"] <- gene_phylotype
 
  gene_haplotype_metrics[gene, "sample_presence"] <- nrow(abun_filt)
  
  gene_haplotype_metrics[gene, "total_num_haplotypes"] <- ncol(abun_filt)

  gene_haplotype_metrics[gene, "min_num_haplotypes"] <- min(num_haplotypes_per_sample)
  gene_haplotype_metrics[gene, "max_num_haplotypes"] <- max(num_haplotypes_per_sample)
  
  gene_haplotype_metrics[gene, "mean_num_haplotypes"] <- mean(num_haplotypes_per_sample)
  gene_haplotype_metrics[gene, "sd_num_haplotypes"] <- sd(num_haplotypes_per_sample)
  gene_haplotype_metrics[gene, "cv_num_haplotypes"] <- sd(num_haplotypes_per_sample) / mean(num_haplotypes_per_sample)
  
  gene_haplotype_metrics[gene, "shannon_mean"] <- mean(gene_shannon_out)
  gene_haplotype_metrics[gene, "shannon_sd"] <- sd(gene_shannon_out)
  gene_haplotype_metrics[gene, "shannon_cv"] <- sd(gene_shannon_out) / mean(gene_shannon_out)
  
  gene_haplotype_metrics[gene, "simpson_mean"] <- mean(gene_simpson_out)
  gene_haplotype_metrics[gene, "simpson_sd"] <- sd(gene_simpson_out)
  gene_haplotype_metrics[gene, "simpson_cv"] <- sd(gene_simpson_out) / mean(gene_simpson_out)
  
  gene_haplotype_metrics[gene, "invsimpson_mean"] <- mean(gene_invsimpson_out)
  gene_haplotype_metrics[gene, "invsimpson_sd"] <- sd(gene_invsimpson_out)
  gene_haplotype_metrics[gene, "invsimpson_cv"] <- sd(gene_invsimpson_out) / mean(gene_invsimpson_out)
  

  
  gene_haplotype_metrics[gene, "mean_norm_num_haplotypes"] <- mean(norm_num_haplotypes_per_sample)
  gene_haplotype_metrics[gene, "sd_norm_num_haplotypes"] <- sd(norm_num_haplotypes_per_sample)
  gene_haplotype_metrics[gene, "cv_norm_num_haplotypes"] <- sd(norm_num_haplotypes_per_sample) / mean(norm_num_haplotypes_per_sample)
  
  gene_haplotype_metrics[gene, "shannon_norm_mean"] <- mean(norm_shannon_per_sample)
  gene_haplotype_metrics[gene, "shannon_norm_sd"] <- sd(norm_shannon_per_sample)
  gene_haplotype_metrics[gene, "shannon_norm_cv"] <- sd(norm_shannon_per_sample) / mean(norm_shannon_per_sample)
  
  gene_haplotype_metrics[gene, "simpson_norm_mean"] <- mean(norm_simpson_per_sample)
  gene_haplotype_metrics[gene, "simpson_norm_sd"] <- sd(norm_simpson_per_sample)
  gene_haplotype_metrics[gene, "simpson_norm_cv"] <- sd(norm_simpson_per_sample) / mean(norm_simpson_per_sample)
  
  gene_haplotype_metrics[gene, "invsimpson_norm_mean"] <- mean(norm_invsimpson_per_sample)
  gene_haplotype_metrics[gene, "invsimpson_norm_sd"] <- sd(norm_invsimpson_per_sample)
  gene_haplotype_metrics[gene, "invsimpson_norm_cv"] <- sd(norm_invsimpson_per_sample) / mean(norm_invsimpson_per_sample)

}

gene_haplotype_metrics$gene <- rownames(gene_haplotype_metrics)

gene_haplotype_metrics$num_ref_genomes <- panaroo_summary$panaroo_all_phylotypes[rownames(gene_haplotype_metrics), "No..isolates"]

gene_haplotype_metrics$percent_ref_genomes <- panaroo_summary$panaroo_all_phylotypes[rownames(gene_haplotype_metrics), "percent_isolates"]

gene_haplotype_metrics <- gene_haplotype_metrics[, c("gene", "num_ref_genomes", "percent_ref_genomes", col_categories)]

write.table(x = gene_haplotype_metrics,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/gene_haplotype_metrics.tsv",
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)
