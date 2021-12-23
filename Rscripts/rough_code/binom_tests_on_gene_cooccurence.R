rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)
library(reshape2)
library(dplyr)
library(ggbeeswarm)
library(plyr)


setwd("honey_bee_pangenome/")
source("Rscripts/functions.R")

path_to_panaroo <- "data/binning_output/panaroo_output_dastools_and_ref/Gilliamella_gene_presence_absence_roary-formatted.csv"

Gilliamella_panaroo_out <- read.table(path_to_panaroo,
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilli", rownames(Gilliamella_panaroo_out), sep = "_")

path_to_basic_breadth <- "data/Ellegaard_panaroo_pangenome_mapping/bowtie2_mapped_gene_coverage/best_hit_simple_comp_map/Gilliamella_coverage_breadth/"

Gilliamella_basic_breadth <- read_in_breadth_files(in_path = path_to_basic_breadth,
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
Gilliamella_present[Gilliamella_present < 0.5] <- 0
Gilliamella_present[Gilliamella_present >= 0.5] <- 1

Gilliamella_genes2test <- Gilliamella_present
Gilliamella_genes2test <- Gilliamella_genes2test[which(rowSums(Gilliamella_genes2test) >= 3), ]
Gilliamella_genes2test <- Gilliamella_genes2test[which(rowSums(Gilliamella_genes2test) <= 0.33 * ncol(Gilliamella_genes2test)), ]
Gilliamella_genes2test <- rownames(Gilliamella_genes2test)

Gilliamella_present_genes2test <- Gilliamella_present[Gilliamella_genes2test, ]

binom_test_p <- c()
binom_test_output <- list()

for (Gilliamella_gene1 in Gilliamella_genes2test) {

  Gilliamella_gene1_pos_samples <- colnames(Gilliamella_present_genes2test)[which(as.numeric(Gilliamella_present_genes2test[Gilliamella_gene1, ]) == 1)]
  Gilliamella_gene1_pos_freq <- length(Gilliamella_gene1_pos_samples) / ncol(Gilliamella_present_genes2test)
  
  
  for (Gilliamella_gene2 in Gilliamella_genes2test) {
  
    if (Gilliamella_gene1 == Gilliamella_gene2) { next }
    
    Gilliamella_gene2_pos_samples <- colnames(Gilliamella_present_genes2test)[which(as.numeric(Gilliamella_present_genes2test[Gilliamella_gene2, ]) == 1)]
    Gilliamella_gene2_pos_freq <- length(Gilliamella_gene2_pos_samples) / ncol(Gilliamella_present_genes2test)
  
    Gilliamella_obs_intersect <- length(which(Gilliamella_gene1_pos_samples %in% Gilliamella_gene2_pos_samples))
    
    Gilliamella_compare_id <- paste(Gilliamella_gene1, Gilliamella_gene2, sep = "|")
    
    binom_test_output[[Gilliamella_compare_id]] <- binom.test(x = Gilliamella_obs_intersect,
                                                              n = length(Gilliamella_gene1_pos_samples),
                                                              p = Gilliamella_gene2_pos_freq)

    binom_test_p <- c(binom_test_p, binom_test_output[[Gilliamella_compare_id]]$p.value)
  
  }
}

bonf_binom_test_p <- p.adjust(binom_test_p, "bonf")

