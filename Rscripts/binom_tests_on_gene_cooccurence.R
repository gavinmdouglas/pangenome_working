### Code to run coincident analysis between all Gilliamella and Snodgrassella genes within a certain window. ###

rm(list = ls(all.names = TRUE))

library("parallel")
library("plyr")

binom_cooccur_test <- function(focal_gene, presence_table) {
  
  binom_test_output <- list()
  
  all_samples <- colnames(presence_table)
  num_samples <- ncol(presence_table)
  
  focal_gene_pos_samples <- all_samples[which(as.numeric(presence_table[focal_gene, ]) == 1)]
  focal_gene_pos_freq <- length(focal_gene_pos_samples) / num_samples
  
  genes2test <- rownames(presence_table)
  genes2test <- genes2test[-c(1:which(genes2test == focal_gene))]
  
  presence_table <- presence_table[genes2test, , drop = FALSE]
  
  for (gene2test in genes2test) {
    
    gene2test_pos_samples <- all_samples[which(as.numeric(presence_table[gene2test, ]) == 1)]
    gene2test_pos_freq <- length(gene2test_pos_samples) / num_samples
    
    obs_intersect <- length(which(focal_gene_pos_samples %in% gene2test_pos_samples))
    
    binom_test_output[[gene2test]] <- binom.test(x = obs_intersect,
                                                 n = num_samples,
                                                 p = focal_gene_pos_freq * gene2test_pos_freq)
  }
  
  return(binom_test_output)
  
}

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


Gilliamella_path_to_panaroo <- "/home/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Gilliamella_gene_presence_absence_roary-formatted.csv.gz"
Gilliamella_path_to_basic_breadth <- "/home/gdouglas/projects/honey_bee/Ellegaard/pangenome_competitive_mapping/bedtools_coverage_out/panaroo_best_hit/Gilliamella_coverage_breadth/"
Snodgrassella_path_to_panaroo <- "/home/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Snodgrassella_gene_presence_absence_roary-formatted.csv.gz"
Snodgrassella_path_to_basic_breadth <- "/home/gdouglas/projects/honey_bee/Ellegaard/pangenome_competitive_mapping/bedtools_coverage_out/panaroo_best_hit/Snodgrassella_coverage_breadth/"
presence_cutoff <- 0.5

Gilliamella_panaroo_out <- read.table(Gilliamella_path_to_panaroo,
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilli", rownames(Gilliamella_panaroo_out), sep = "_")


Gilliamella_basic_breadth <- read_in_breadth_files(in_path = Gilliamella_path_to_basic_breadth,
                                                     pattern = ".breadth.bedGraph.gz",
                                                     roary_formatted_pangenome = Gilliamella_panaroo_out)

Gilliamella_breadth_by_sample_clean <- list()

for (Gilliamella_sample in names(Gilliamella_basic_breadth$breadth_by_sample)) {
  Gilliamella_breadth_by_sample_clean[[Gilliamella_sample]] <- Gilliamella_basic_breadth$breadth_by_sample[[Gilliamella_sample]][, c("gene", "V7"), drop = FALSE]
  colnames(Gilliamella_breadth_by_sample_clean[[Gilliamella_sample]]) <- c("gene", Gilliamella_sample)
}

Gilliamella_present <- join_all(Gilliamella_breadth_by_sample_clean, by="gene", type='left')
rownames(Gilliamella_present) <- Gilliamella_present$gene
Gilliamella_present <- Gilliamella_present[, -which(colnames(Gilliamella_present) == "gene")]
Gilliamella_present[Gilliamella_present < presence_cutoff] <- 0
Gilliamella_present[Gilliamella_present >= presence_cutoff] <- 1

Gilliamella_noncore_genes <- rownames(Gilliamella_panaroo_out)[which(Gilliamella_panaroo_out$No..isolates < 0.95 * max(Gilliamella_panaroo_out$No..isolates))]

Gilliamella_present_noncore <- Gilliamella_present[which(rownames(Gilliamella_present) %in% Gilliamella_noncore_genes), ]




Snodgrassella_panaroo_out <- read.table(Snodgrassella_path_to_panaroo,
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snod", rownames(Snodgrassella_panaroo_out), sep = "_")

Snodgrassella_basic_breadth <- read_in_breadth_files(in_path = Snodgrassella_path_to_basic_breadth,
                                                   pattern = ".breadth.bedGraph.gz",
                                                   roary_formatted_pangenome = Snodgrassella_panaroo_out)

Snodgrassella_breadth_by_sample_clean <- list()

for (Snodgrassella_sample in names(Snodgrassella_basic_breadth$breadth_by_sample)) {
  Snodgrassella_breadth_by_sample_clean[[Snodgrassella_sample]] <- Snodgrassella_basic_breadth$breadth_by_sample[[Snodgrassella_sample]][, c("gene", "V7"), drop = FALSE]
  colnames(Snodgrassella_breadth_by_sample_clean[[Snodgrassella_sample]]) <- c("gene", Snodgrassella_sample)
}

Snodgrassella_present <- join_all(Snodgrassella_breadth_by_sample_clean, by="gene", type='left')
rownames(Snodgrassella_present) <- Snodgrassella_present$gene
Snodgrassella_present <- Snodgrassella_present[, -which(colnames(Snodgrassella_present) == "gene")]
Snodgrassella_present[Snodgrassella_present < presence_cutoff] <- 0
Snodgrassella_present[Snodgrassella_present >= presence_cutoff] <- 1

Snodgrassella_noncore_genes <- rownames(Snodgrassella_panaroo_out)[which(Snodgrassella_panaroo_out$No..isolates < 0.95 * max(Snodgrassella_panaroo_out$No..isolates))]

Snodgrassella_present_noncore <- Snodgrassella_present[which(rownames(Snodgrassella_present) %in% Snodgrassella_noncore_genes), ]

Gilliamella_and_Snodgrassella_present_noncore <- rbind(Gilliamella_present_noncore, Snodgrassella_present_noncore)

Gilliamella_and_Snodgrassella_present_noncore <- Gilliamella_and_Snodgrassella_present_noncore[which(rowSums(Gilliamella_and_Snodgrassella_present_noncore) >= 5), ]
Gilliamella_and_Snodgrassella_present_noncore <- Gilliamella_and_Snodgrassella_present_noncore[which(rowSums(Gilliamella_and_Snodgrassella_present_noncore) <= 0.75 * ncol(Gilliamella_and_Snodgrassella_present_noncore)), ]

# Remove unneeded intermediate files.
rm(Gilliamella_basic_breadth)
rm(Gilliamella_breadth_by_sample_clean)
rm(Gilliamella_present)
rm(Snodgrassella_basic_breadth)
rm(Snodgrassella_breadth_by_sample_clean)
rm(Snodgrassella_present)

binom_test_output <- mclapply(rownames(Gilliamella_and_Snodgrassella_present_noncore),
                              FUN = binom_cooccur_test,
                              presence_table = Gilliamella_and_Snodgrassella_present_noncore,
                              mc.cores = 25)

names(binom_test_output) <- rownames(Gilliamella_and_Snodgrassella_present_noncore)

saveRDS(object = binom_test_output,
        file = "/home/gdouglas/projects/honey_bee/Ellegaard/coincidence_analysis/Gilliamella_and_Snodgrassella_panaroo_coincidence_binomial.rds")
