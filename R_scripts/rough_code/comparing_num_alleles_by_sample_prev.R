rm(list = ls(all.names = TRUE))

library(extraDistr)

all_genes_abun <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/all_genes_OTU.rds")

num_alleles_per_gene <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/num_alleles_per_gene.tsv",
                              header = FALSE, sep = "\t", row.names = 1)
Bifidobacterium_asteroides_num_alleles <- num_alleles_per_gene[grep("Bifidobacterium_asteroides", rownames(num_alleles_per_gene)), , drop = FALSE]
Bifidobacterium_asteroides_invariant <- rownames(Bifidobacterium_asteroides_num_alleles)[which(Bifidobacterium_asteroides_num_alleles$V2 == 1)]
hist(all_genes_abun$Bartonella_apis$Bartonella_apis_aaeA)

hist(rowSums(all_genes_abun$Bifidobacterium_asteroides$`Bifidobacterium_asteroides_lacZ~~~lacZ_2~~~lacZ_1` > 0))

test <- rowSums(all_genes_abun$Bifidobacterium_asteroides$`Bifidobacterium_asteroides_lacZ~~~lacZ_2~~~lacZ_1` > 0)
data_pl <- displ$new(test)
est <- estimate_xmin(data_pl)
data_pl$xmin <- est$xmin
data_pl$pars <- est$pars

sample_prev <- sapply(all_genes_abun$Bifidobacterium_asteroides, nrow)

min_allele_count <- sapply(all_genes_abun$Bifidobacterium_asteroides, function(x) { min(rowSums(x > 0)) })

mean_allele_count <- sapply(all_genes_abun$Bifidobacterium_asteroides, function(x) { mean(rowSums(x > 0)) })

Bifidobacterium_asteroides_invariant_missing <- Bifidobacterium_asteroides_invariant[which(! Bifidobacterium_asteroides_invariant %in% names(all_genes_abun$Bifidobacterium_asteroides))]

missing_sample_prev <- c()
for (g in Bifidobacterium_asteroides_invariant_missing) {
  sample_file <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/starting_files/", g, "_samples.txt", sep = "")
  missing_sample_prev <- c(missing_sample_prev, nrow(read.table(sample_file)))
}

prev_by_allele_count <- data.frame(gene = names(all_genes_abun$Bifidobacterium_asteroides),
                                   prevalence = sample_prev,
                                   min_allele_count = min_allele_count,
                                   mean_allele_count = mean_allele_count)

prev_by_allele_count_invaraint_missing <- data.frame(gene = Bifidobacterium_asteroides_invariant_missing,
                                                     prevalence = missing_sample_prev,
                                                     min_allele_count = 1,
                                                     mean_allele_count = 1)

prev_by_allele_count <- rbind(prev_by_allele_count, prev_by_allele_count_invaraint_missing)

max_prev <- max(prev_by_allele_count$prevalence)
max_min_allele <- max(prev_by_allele_count$min_allele_count)
num_comparisons <- max_prev * max_min_allele

prev_by_allele_count_for_heatmap <- data.frame(matrix(NA, nrow = num_comparisons, ncol = 3))
colnames(prev_by_allele_count_for_heatmap) <- c("prevalence", "min_allele_count", "count")

comparison_number <- 1
for (i in 1:max_prev) {
  
  for (j in 1:max_min_allele) {
    
    prev_by_allele_count_for_heatmap[comparison_number, ] <- c(i, j, length(which(prev_by_allele_count$prevalence == i & prev_by_allele_count$min_allele_count == j)))
    
    comparison_number <- comparison_number + 1
  }
}

library(ggplot2)

prev_by_allele_count_for_heatmap$count[which(prev_by_allele_count_for_heatmap$count == 0)] <- NA

ggplot(data = prev_by_allele_count_for_heatmap, aes(x = prevalence, y = min_allele_count, fill = log10(count))) +
  geom_tile() +
  theme_bw()

obs_counts <- rowSums(all_genes_abun$Bifidobacterium_asteroides$`Bifidobacterium_asteroides_lacZ~~~lacZ_2~~~lacZ_1` > 0)

cleaned_obs_counts <- c()
exp_prob <- c()

obs_mean <- mean(obs_counts)

for (i in 1:50) {
  cleaned_obs_counts <- c(cleaned_obs_counts, length(which(obs_counts == i)))
  exp_prob <- c(exp_prob, dtpois(i, lambda = obs_mean))
}

exp_prob <- exp_prob / sum(exp_prob)

ks.test(x = cleaned_obs_counts, y = exp_prob)

chisq.test(x = cleaned_obs_counts, p = exp_prob)
