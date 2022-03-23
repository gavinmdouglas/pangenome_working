identify_enriched_categories <- function(genes,
                                         background,
                                         category_to_gene_map) {
  
  enrichments_out <- data.frame(matrix(NA, nrow = length(category_to_gene_map), ncol = 8))
  rownames(enrichments_out) <- names(category_to_gene_map)
  colnames(enrichments_out) <- c("category", "genes_num_category", "genes_num_other",
                                 "background_num_category", "background_num_other", "OR", "p", "fdr")
  
  enrichments_out[names(category_to_gene_map), "category"] <- names(category_to_gene_map)
  
  for (category in rownames(enrichments_out)) {
    
    genes_num_category <- length(which(genes %in% category_to_gene_map[[category]]))
    genes_num_other <- length(genes) - genes_num_category
    
    background_num_category <- length(which(background %in% category_to_gene_map[[category]]))
    background_num_other <- length(background) - background_num_category
    
    count_table <- matrix(c(genes_num_category, genes_num_other, background_num_category, background_num_other), nrow = 2, ncol = 2)
    
    fisher_out <- fisher.test(count_table)
    
    enrichments_out[category, c("genes_num_category",
                                "genes_num_other",
                                "background_num_category",
                                "background_num_other", "p")] <- c(genes_num_category,
                                                                   genes_num_other,
                                                                   background_num_category,
                                                                   background_num_other,
                                                                   fisher_out$p.value)
    if (genes_num_other > 0) {
      ratio_numer <- genes_num_category / genes_num_other
    } else {
      ratio_numer <- genes_num_category / 1 
    }
    
    if (background_num_other == 0) {
      ratio_denom <- 1
    } else if(background_num_category == 0) {
      ratio_denom <- 1 / background_num_other
    } else {
      ratio_denom <- background_num_category / background_num_other
    }
    
    enrichments_out[category, "OR"] <- ratio_numer / ratio_denom
  }
  
  enrichments_out$fdr <- p.adjust(enrichments_out$p, "fdr")
  
  rownames(enrichments_out) <- NULL
  
  return(enrichments_out)
  
}



read_in_breadth_files <- function(in_path, pattern, roary_formatted_pangenome) {
  
  in_path <- gsub("/$", "", in_path)
  
  # Read in all breadth files.
  input_breadth_files <- list.files(path = in_path, full.names = TRUE, pattern = pattern)
  
  input_breadth_by_sample <- lapply(input_breadth_files, function(x) { read.table(x, header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE) })
  input_samples <- gsub(paste(in_path, "/", sep = ""), "", input_breadth_files)
  #input_samples <- gsub("\\..*\\.merged.nonparalog.breadth.bedGraph.gz", "", input_samples)
  input_samples <- gsub("\\.breadth.bedGraph.gz", "", input_samples)
  names(input_breadth_by_sample) <- input_samples
  
  input_breadth_by_sample[[1]]$sample <- input_samples[1]
  input_breadth_by_sample[[1]]$gene <- rownames(input_breadth_by_sample[[1]])
  rownames(input_breadth_by_sample[[1]]) <- NULL
  input_breadth_combined <- input_breadth_by_sample[[1]]
  
  for (i in 2:length(input_samples)) {
    input_breadth_by_sample[[i]]$sample <- input_samples[i]
    input_breadth_by_sample[[i]]$gene <- rownames(input_breadth_by_sample[[i]])
    rownames(input_breadth_by_sample[[i]]) <- NULL
    input_breadth_combined <- rbind(input_breadth_combined,  input_breadth_by_sample[[i]])
  }
  
  input_breadth_combined <- input_breadth_combined[, c("sample", "gene", "V7")]
  colnames(input_breadth_combined)[3] <- "breadth"
  
  input_unique_genes <- input_breadth_combined$gene[-which(duplicated(input_breadth_combined$gene))]
  
  # Get per-gene summary across samples.
  input_breadth_summary <- data.frame(matrix(NA,
                                             nrow = length(input_unique_genes),
                                             ncol = 13))
  colnames(input_breadth_summary) <- c("ref_num", "mean", "median", "max", "min", "sd", "cv", "num_at_least_0.1",
                                       "num_at_least_0.25", "num_at_least_0.5", "num_at_least_0.9", "num_1", "num_0")
  rownames(input_breadth_summary) <- input_unique_genes
  
  for (input_gene in input_unique_genes) {
    input_breadth_summary[input_gene, "ref_num"] <- roary_formatted_pangenome[input_gene, "No..isolates"]
    input_gene_breadth_values <- input_breadth_combined[which(input_breadth_combined$gene == input_gene), "breadth"]
    
    input_breadth_summary[input_gene, "mean"] <- mean(input_gene_breadth_values)
    input_breadth_summary[input_gene, "median"] <- median(input_gene_breadth_values)
    input_breadth_summary[input_gene, "max"] <- max(input_gene_breadth_values)
    input_breadth_summary[input_gene, "min"] <- min(input_gene_breadth_values)
    input_breadth_summary[input_gene, "sd"] <- sd(input_gene_breadth_values)
    input_breadth_summary[input_gene, "cv"] <- input_breadth_summary[input_gene, "sd"] / input_breadth_summary[input_gene, "mean"]
    input_breadth_summary[input_gene, "num_at_least_0.1"] <- length(which(input_gene_breadth_values >= 0.1))
    input_breadth_summary[input_gene, "num_at_least_0.25"] <- length(which(input_gene_breadth_values >= 0.25))
    input_breadth_summary[input_gene, "num_at_least_0.5"] <- length(which(input_gene_breadth_values >= 0.5))
    input_breadth_summary[input_gene, "num_at_least_0.9"] <- length(which(input_gene_breadth_values >= 0.9))
    input_breadth_summary[input_gene, "num_1"] <- length(which(input_gene_breadth_values == 1))
    input_breadth_summary[input_gene, "num_0"] <- length(which(input_gene_breadth_values == 0))
  }
  
  # Create rounded version of metrics.
  input_breadth_summary$mean_rounded <- round(input_breadth_summary$mean, digits = 1)
  input_breadth_summary$median_rounded <- round(input_breadth_summary$median, digits = 1)
  input_breadth_summary$max_rounded <- round(input_breadth_summary$max, digits = 1)
  input_breadth_summary$min_rounded <- round(input_breadth_summary$min, digits = 1)
  input_breadth_summary$sd_rounded <- round(input_breadth_summary$sd, digits = 1)
  input_breadth_summary$cv_rounded <- round(input_breadth_summary$cv, digits = 1)
  
  return(list(breadth_by_sample = input_breadth_by_sample,
              breadth_all_samples = input_breadth_combined,
              breadth_summary = input_breadth_summary,
              unique_genes = input_unique_genes))
  
}



read_in_depth_files <- function(in_path, pattern, roary_formatted_pangenome) {
  
  in_path <- gsub("/$", "", in_path)
  
  # Read in all depth files.
  input_depth_files <- list.files(path = in_path, full.names = TRUE, pattern = pattern)
  
  input_depth_by_sample <- lapply(input_depth_files, function(x) { read.table(x, header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE) })
  input_samples <- gsub(paste(in_path, "/", sep = ""), "", input_depth_files)
  input_samples <- gsub("\\.mean.bedGraph.gz", "", input_samples)
  names(input_depth_by_sample) <- input_samples
  
  input_depth_by_sample[[1]]$sample <- input_samples[1]
  input_depth_by_sample[[1]]$gene <- rownames(input_depth_by_sample[[1]])
  rownames(input_depth_by_sample[[1]]) <- NULL
  input_depth_combined <- input_depth_by_sample[[1]]
  
  for (i in 2:length(input_samples)) {
    input_depth_by_sample[[i]]$sample <- input_samples[i]
    input_depth_by_sample[[i]]$gene <- rownames(input_depth_by_sample[[i]])
    rownames(input_depth_by_sample[[i]]) <- NULL
    input_depth_combined <- rbind(input_depth_combined,  input_depth_by_sample[[i]])
  }
  
  input_depth_combined <- input_depth_combined[, c("sample", "gene", "V4")]
  colnames(input_depth_combined)[3] <- "depth"
  
  input_unique_genes <- input_depth_combined$gene[-which(duplicated(input_depth_combined$gene))]
  
  # Get per-gene summary across samples.
  input_depth_summary <- data.frame(matrix(NA,
                                           nrow = length(input_unique_genes),
                                           ncol = 13))
  colnames(input_depth_summary) <- c("ref_num", "mean", "median", "max", "min", "sd", "cv", "num_at_least_0.1",
                                     "num_at_least_0.25", "num_at_least_0.5", "num_at_least_0.9", "num_1", "num_0")
  rownames(input_depth_summary) <- input_unique_genes
  
  for (input_gene in input_unique_genes) {
    input_depth_summary[input_gene, "ref_num"] <- roary_formatted_pangenome[input_gene, "No..isolates"]
    input_gene_depth_values <- input_depth_combined[which(input_depth_combined$gene == input_gene), "depth"]
    
    input_depth_summary[input_gene, "mean"] <- mean(input_gene_depth_values)
    input_depth_summary[input_gene, "median"] <- median(input_gene_depth_values)
    input_depth_summary[input_gene, "max"] <- max(input_gene_depth_values)
    input_depth_summary[input_gene, "min"] <- min(input_gene_depth_values)
    input_depth_summary[input_gene, "sd"] <- sd(input_gene_depth_values)
    input_depth_summary[input_gene, "cv"] <- input_depth_summary[input_gene, "sd"] / input_depth_summary[input_gene, "mean"]
    input_depth_summary[input_gene, "num_at_least_0.1"] <- length(which(input_gene_depth_values >= 0.1))
    input_depth_summary[input_gene, "num_at_least_0.25"] <- length(which(input_gene_depth_values >= 0.25))
    input_depth_summary[input_gene, "num_at_least_0.5"] <- length(which(input_gene_depth_values >= 0.5))
    input_depth_summary[input_gene, "num_at_least_0.9"] <- length(which(input_gene_depth_values >= 0.9))
    input_depth_summary[input_gene, "num_1"] <- length(which(input_gene_depth_values == 1))
    input_depth_summary[input_gene, "num_0"] <- length(which(input_gene_depth_values == 0))
  }
  
  # Create rounded version of metrics.
  input_depth_summary$mean_rounded <- round(input_depth_summary$mean, digits = 1)
  input_depth_summary$median_rounded <- round(input_depth_summary$median, digits = 1)
  input_depth_summary$max_rounded <- round(input_depth_summary$max, digits = 1)
  input_depth_summary$min_rounded <- round(input_depth_summary$min, digits = 1)
  input_depth_summary$sd_rounded <- round(input_depth_summary$sd, digits = 1)
  input_depth_summary$cv_rounded <- round(input_depth_summary$cv, digits = 1)
  
  return(list(depth_by_sample = input_depth_by_sample,
              depth_all_samples = input_depth_combined,
              depth_summary = input_depth_summary,
              unique_genes = input_unique_genes))
  
}



read_KO_pathway_map <- function(filename, min_num_funcs = 1) {
  
  KO_pathway_map <- read.table(filename, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  KO_pathway_map$V1 <- gsub("ko:", "", KO_pathway_map$V1)
  
  KO_pathway_map <- KO_pathway_map[-grep("path:map", KO_pathway_map$V2), ]
  
  KO_pathway_map$V2 <- gsub("path:", "", KO_pathway_map$V2)
  
  pathway_to_KO <- list()
  
  all_pathways <- KO_pathway_map$V2
  all_pathways <- all_pathways[-which(duplicated(all_pathways))]
  
  for (pathway in all_pathways) {
    pathway_funcs <- sort(unique(KO_pathway_map[which(KO_pathway_map$V2 == pathway), "V1"]))
    
    if (length(pathway_funcs) > min_num_funcs) {
      pathway_to_KO[[pathway]] <- pathway_funcs
    }
    
  }
  
  return(pathway_to_KO)
}


read_KO_module_map <- function(filename, min_num_funcs=1) {
  
  KO_module_map <- read.table(filename, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  KO_module_map$V1 <- gsub("ko:", "", KO_module_map$V1)
  
  KO_module_map$V2 <- gsub("md:", "", KO_module_map$V2)
  
  module_to_KO <- list()
  
  all_modules <- KO_module_map$V2
  all_modules <- all_modules[-which(duplicated(all_modules))]
  
  for (module in all_modules) {
    module_funcs <- sort(unique(KO_module_map[which(KO_module_map$V2 == module), "V1"]))
    
    if (length(module_funcs) > min_num_funcs) {
      module_to_KO[[module]] <- module_funcs
    }
  }
  
  return(module_to_KO)
}



test_for_enriched_gene_families <- function(gene_set, background_set, min_count_in_background = 1) {
  
  gene_set_multi_i <- grep(",", gene_set)
  if (length(gene_set_multi_i) > 0) {
    
    gene_set_multi <- gene_set[gene_set_multi_i]
    gene_set <- gene_set[-gene_set_multi_i]
    
    for (multi_gene in gene_set_multi) {
      gene_set <- c(gene_set, str_split(multi_gene, ",")[[1]])
    }
    
  }
  
  background_set_multi_i <- grep(",", background_set)
  if (length(background_set_multi_i) > 0) {
    
    background_set_multi <- background_set[background_set_multi_i]
    background_set <- background_set[-background_set_multi_i]
    
    for (multi_background_gene in background_set_multi) {
      background_set <- c(background_set, str_split(multi_background_gene, ",")[[1]])
    }
    
  }
  
  background_set_breakdown <- table(background_set)
  
  gene_families_to_test <- names(background_set_breakdown)[which(background_set_breakdown >= min_count_in_background)]
  
  if (length(gene_families_to_test) == 0) { return(NA) }
  
  fisher_test_p <- c()
  fisher_test_or <- c()
  
  for (gene_family in gene_families_to_test) {
    gene_set_gene_family_count <- length(which(gene_set == gene_family))
    gene_set_non_gene_family_count <- length(which(gene_set != gene_family)) 
    background_set_family_count <- length(which(background_set == gene_family)) 
    background_set_non_family_count <- length(which(background_set != gene_family)) 
    
    fisher_test_output <- fisher.test(matrix(c(gene_set_gene_family_count, gene_set_non_gene_family_count,
                                               background_set_family_count, background_set_non_family_count),
                                             ncol = 2))
    
    fisher_test_p <- c(fisher_test_p, fisher_test_output$p.value)
    fisher_test_or <- c(fisher_test_or, fisher_test_output$estimate)
    
  }
  
  names(fisher_test_p) <- gene_families_to_test
  names(fisher_test_or) <- gene_families_to_test
  
  fisher_test_fdr <- p.adjust(fisher_test_p, "BH")
  
  return(list(fisher_p=fisher_test_p,
              fisher_BH=fisher_test_fdr,
              fisher_or=fisher_test_or))
  
}

