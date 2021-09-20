read_in_breadth_files <- function(in_path, pattern, roary_formatted_pangenome) {
  
  in_path <- gsub("/$", "", in_path)
  
  # Read in all breadth files.
  input_breadth_files <- list.files(path = in_path, full.names = TRUE, pattern = pattern)
  
  input_breadth_by_sample <- lapply(input_breadth_files, function(x) { read.table(x, header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE) })
  input_samples <- gsub(paste(in_path, "/", sep = ""), "", input_breadth_files)
  input_samples <- gsub("\\..*\\.merged.nonparalog.breadth.bedGraph.gz", "", input_samples)
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
  input_samples <- gsub("\\..*\\.merged.nonparalog.mean.bedGraph.gz", "", input_samples)
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