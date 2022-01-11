rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)
library(reshape2)
library(dplyr)
library(ggbeeswarm)

setwd("honey_bee_pangenome/data/Ellegaard_panaroo_pangenome_mapping/bowtie2_mapped_gene_coverage/best_hit_simple_comp_map/")

# Read in reference genome info.
Snodgrassella_panaroo_out <- read.table("../../../binning_output/panaroo_output_dastools_and_ref/Snodgrassella_gene_presence_absence_roary-formatted.csv",
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snod", rownames(Snodgrassella_panaroo_out), sep = "_")

# Read in all breadth files.
Snodgrassella_breadth_files <- list.files(path = "Snodgrassella_coverage_breadth/",
                                        full.names = TRUE, pattern = ".breadth.bedGraph.gz")

Snodgrassella_breadth_by_sample <- lapply(Snodgrassella_breadth_files, function(x) { read.table(x, header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE) })
Snodgrassella_samples <- gsub("Snodgrassella_coverage_breadth//", "", Snodgrassella_breadth_files)
Snodgrassella_samples <- gsub(".Snodgrassella.merged.nonparalog.breadth.bedGraph.gz", "", Snodgrassella_samples)
names(Snodgrassella_breadth_by_sample) <- Snodgrassella_samples

# Quick breadth histograms.
par(mfrow = c(1, 2))
hist(Snodgrassella_breadth_by_sample[[13]]$V7, breaks = 100)
hist(Snodgrassella_breadth_by_sample[[32]]$V7, breaks = 100)
par(mfrow = c(1, 1))

Snodgrassella_breadth_by_sample[[1]]$sample <- Snodgrassella_samples[1]
Snodgrassella_breadth_by_sample[[1]]$gene <- rownames(Snodgrassella_breadth_by_sample[[1]])
rownames(Snodgrassella_breadth_by_sample[[1]]) <- NULL
Snodgrassella_breadth_combined <- Snodgrassella_breadth_by_sample[[1]]
for (i in 2:length(Snodgrassella_samples)) {
  Snodgrassella_breadth_by_sample[[i]]$sample <- Snodgrassella_samples[i]
  Snodgrassella_breadth_by_sample[[i]]$gene <- rownames(Snodgrassella_breadth_by_sample[[i]])
  rownames(Snodgrassella_breadth_by_sample[[i]]) <- NULL
  Snodgrassella_breadth_combined <- rbind(Snodgrassella_breadth_combined,  Snodgrassella_breadth_by_sample[[i]])
}

Snodgrassella_breadth_combined <- Snodgrassella_breadth_combined[, c("sample", "gene", "V7")]
colnames(Snodgrassella_breadth_combined)[3] <- "breadth"

# Combined breadth histogram across all samples (i.e., each gene represented 54 times)
hist(Snodgrassella_breadth_combined$breadth)
# The above shows that whether to call genes with any coverage at all as present or not could have a big impact!

Snodgrassella_unique_genes <- Snodgrassella_breadth_combined$gene[-which(duplicated(Snodgrassella_breadth_combined$gene))]

# Get per-gene summary across samples.
Snodgrassella_breadth_summary <- data.frame(matrix(NA,
                                                 nrow = length(Snodgrassella_unique_genes),
                                                 ncol = 12))
colnames(Snodgrassella_breadth_summary) <- c("ref_num", "mean", "median", "max", "sd", "cv", "num_at_least_0.1",
                                           "num_at_least_0.25", "num_at_least_0.5", "num_at_least_0.9", "num_1", "num_0")
rownames(Snodgrassella_breadth_summary) <- Snodgrassella_unique_genes

for (Snodgrassella_gene in Snodgrassella_unique_genes) {
  Snodgrassella_breadth_summary[Snodgrassella_gene, "ref_num"] <- Snodgrassella_panaroo_out[Snodgrassella_gene, "No..isolates"]
  Snodgrassella_gene_breadth_values <- Snodgrassella_breadth_combined[which(Snodgrassella_breadth_combined$gene == Snodgrassella_gene), "breadth"]
  
  Snodgrassella_breadth_summary[Snodgrassella_gene, "mean"] <- mean(Snodgrassella_gene_breadth_values)
  Snodgrassella_breadth_summary[Snodgrassella_gene, "median"] <- median(Snodgrassella_gene_breadth_values)
  Snodgrassella_breadth_summary[Snodgrassella_gene, "max"] <- max(Snodgrassella_gene_breadth_values)
  Snodgrassella_breadth_summary[Snodgrassella_gene, "sd"] <- sd(Snodgrassella_gene_breadth_values)
  Snodgrassella_breadth_summary[Snodgrassella_gene, "cv"] <- Snodgrassella_breadth_summary[Snodgrassella_gene, "sd"] / Snodgrassella_breadth_summary[Snodgrassella_gene, "mean"]
  Snodgrassella_breadth_summary[Snodgrassella_gene, "num_at_least_0.1"] <- length(which(Snodgrassella_gene_breadth_values >= 0.1))
  Snodgrassella_breadth_summary[Snodgrassella_gene, "num_at_least_0.25"] <- length(which(Snodgrassella_gene_breadth_values >= 0.25))
  Snodgrassella_breadth_summary[Snodgrassella_gene, "num_at_least_0.5"] <- length(which(Snodgrassella_gene_breadth_values >= 0.5))
  Snodgrassella_breadth_summary[Snodgrassella_gene, "num_at_least_0.9"] <- length(which(Snodgrassella_gene_breadth_values >= 0.9))
  Snodgrassella_breadth_summary[Snodgrassella_gene, "num_1"] <- length(which(Snodgrassella_gene_breadth_values == 1))
  Snodgrassella_breadth_summary[Snodgrassella_gene, "num_0"] <- length(which(Snodgrassella_gene_breadth_values == 0))
}


# Mean breadth vs num. ref genomes with gene
Snodgrassella_breadth_summary$mean_rounded <- round(Snodgrassella_breadth_summary$mean, digits = 1)

group_by(Snodgrassella_breadth_summary, ref_num, mean_rounded) %>%
  summarize(n = log(n())) %>% 
  ggplot(aes(ref_num, mean_rounded, fill = n)) +
  geom_tile() +
  xlab("No. ref. genomes encoding") +
  ylab("Mean breadth of coverage (rounded)") +
  labs(fill='log10(Count)')

# Can do the same for median, sd, etc. etc.

# Heatmaps of no. samples w gene above breadth cut-off vs num ref genomes that encode gene
group_by(Snodgrassella_breadth_summary, ref_num, num_at_least_0.1) %>%
          summarize(n = log(n())) %>% 
          ggplot(aes(ref_num, num_at_least_0.1, fill = n)) +
          geom_tile() +
          xlab("No. ref. genomes encoding") +
          ylab("No. samples with coverage >= 10%") +
          labs(fill='log10(Count)')

group_by(Snodgrassella_breadth_summary, ref_num, num_at_least_0.9) %>%
        summarize(n = log(n())) %>% 
        ggplot(aes(ref_num, num_at_least_0.9, fill = n)) +
        geom_tile() +
        xlab("No. ref. genomes encoding") +
        ylab("No. samples with coverage >= 90%") +
        labs(fill='log10(Count)')


# Exclude genes where coverage is never above 10%
# These could be impossible to map to unambiguously (or could simply not be present at all, but hard to say).
Snodgrassella_negligible_genes <- rownames(Snodgrassella_breadth_summary)[which(Snodgrassella_breadth_summary$num_at_least_0.1 == 0)]
Snodgrassella_nonnegligible_genes <- rownames(Snodgrassella_breadth_summary)[which(Snodgrassella_breadth_summary$num_at_least_0.1 > 0)]

Snodgrassella_breadth_summary_nonnegligible <- Snodgrassella_breadth_summary[Snodgrassella_nonnegligible_genes, ]
Snodgrassella_panaroo_out_nonnegligible_subset <- Snodgrassella_panaroo_out[Snodgrassella_nonnegligible_genes, ]

Snodgrassella_panaroo_nonnegligible_soft_core <- rownames(Snodgrassella_panaroo_out_nonnegligible_subset)[which(Snodgrassella_panaroo_out_nonnegligible_subset$No..isolates >= 0.95 * max(Snodgrassella_panaroo_out_nonnegligible_subset$No..isolates))]

# Soft core genes have lower CV compared to other non-negligible genes.
Snodgrassella_panaroo_nonnegligible_NONsoft_core <- rownames(Snodgrassella_breadth_summary_nonnegligible)[which(! rownames(Snodgrassella_breadth_summary_nonnegligible) %in% Snodgrassella_panaroo_nonnegligible_soft_core)]
boxplot(Snodgrassella_breadth_summary_nonnegligible[Snodgrassella_panaroo_nonnegligible_soft_core, "cv"],
        Snodgrassella_breadth_summary_nonnegligible[Snodgrassella_panaroo_nonnegligible_NONsoft_core, "cv"])

# Table of mean + SD percent soft core genes (>=95% ref. genomes) at 10%, 25%, 50%, and 90% coverage
Snodgrassella_soft_core_coverage_summary_tab <- data.frame(matrix(NA, nrow = 4, ncol = 2))
rownames(Snodgrassella_soft_core_coverage_summary_tab) <- c(">=10%", ">=25%", ">=50%", ">=90%")
colnames(Snodgrassella_soft_core_coverage_summary_tab) <- c("mean", "sd")

Snodgrassella_count_10per <- c()
Snodgrassella_count_25per <- c()
Snodgrassella_count_50per <- c()
Snodgrassella_count_90per <- c()

for (Snodgrassella_sample in Snodgrassella_samples) {
  Snodgrassella_tmp_sample_soft_core_breadth <- Snodgrassella_breadth_by_sample[[Snodgrassella_sample]][which(Snodgrassella_breadth_by_sample[[Snodgrassella_sample]]$gene %in% Snodgrassella_panaroo_nonnegligible_soft_core), ]
  Snodgrassella_count_10per <- c(Snodgrassella_count_10per, (length(which(Snodgrassella_tmp_sample_soft_core_breadth$V7 >= 0.1)) / length(Snodgrassella_panaroo_nonnegligible_soft_core)) * 100)
  Snodgrassella_count_25per <- c(Snodgrassella_count_25per, (length(which(Snodgrassella_tmp_sample_soft_core_breadth$V7 >= 0.25)) / length(Snodgrassella_panaroo_nonnegligible_soft_core)) * 100)
  Snodgrassella_count_50per <- c(Snodgrassella_count_50per, (length(which(Snodgrassella_tmp_sample_soft_core_breadth$V7 >= 0.5)) / length(Snodgrassella_panaroo_nonnegligible_soft_core)) * 100)
  Snodgrassella_count_90per <- c(Snodgrassella_count_90per, (length(which(Snodgrassella_tmp_sample_soft_core_breadth$V7 >= 0.9)) / length(Snodgrassella_panaroo_nonnegligible_soft_core)) * 100)
}

Snodgrassella_soft_core_coverage_summary_tab[">=10%", ] <- c(mean(Snodgrassella_count_10per), sd(Snodgrassella_count_10per))
Snodgrassella_soft_core_coverage_summary_tab[">=25%", ] <- c(mean(Snodgrassella_count_25per), sd(Snodgrassella_count_25per))
Snodgrassella_soft_core_coverage_summary_tab[">=50%", ] <- c(mean(Snodgrassella_count_50per), sd(Snodgrassella_count_50per))
Snodgrassella_soft_core_coverage_summary_tab[">=90%", ] <- c(mean(Snodgrassella_count_90per), sd(Snodgrassella_count_90per))

Snodgrassella_soft_core_percent_coverage <- data.frame(sample = Snodgrassella_samples,
                                                     per10 = Snodgrassella_count_10per,
                                                     per25 = Snodgrassella_count_25per,
                                                     per50 = Snodgrassella_count_50per,
                                                     per90 = Snodgrassella_count_90per)

Snodgrassella_soft_core_percent_coverage_melt <- melt(Snodgrassella_soft_core_percent_coverage, id = "sample")

set.seed(131)
ggplot(data = Snodgrassella_soft_core_percent_coverage_melt, aes(x = variable, y = value)) +
       geom_beeswarm(cex = 2) +
       geom_boxplot(outlier.shape = NA, fill = "grey", alpha = 0.4) +
       theme_bw() +
       ylab("% of (non-negligible) soft core genes called as present") +
       xlab("% breadth cut-off used to call a gene present")


# Exploring how binomial test would work (with binom.test)

test_data = matrix(sample(c(0, 1), replace=TRUE, size = 20), nrow = 2, ncol = 10)

num_rep <- 100000
num_success = 0
for (rep in 1:num_rep) {
  # Gene to randomly observe (first in vector - unobserved is second)
  ran_gene <- sample(c(1, 2), size = 2)
  
  observed_genome_index <- sample(which(test_data[ran_gene[1], ] == 1), 1)
  
  num_success <- num_success + test_data[ran_gene[2], observed_genome_index]
}

num_success / num_rep
