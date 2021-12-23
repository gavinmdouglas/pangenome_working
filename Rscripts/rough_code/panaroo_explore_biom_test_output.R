rm(list = ls(all.names = TRUE))

panaroo_coincidence_binomial_raw <- readRDS("/home/gdouglas/projects/honey_bee/Ellegaard/coincidence_analysis/Gilliamella_and_Snodgrassella_panaroo_coincidence_binomial.rds")

total_num_comparisons <- 0

for (i in 1:length(panaroo_coincidence_binomial_raw)) {
  total_num_comparisons <- total_num_comparisons + length(panaroo_coincidence_binomial_raw[[i]])
}

panaroo_coincidence_binomial <- data.frame(matrix(NA, nrow = total_num_comparisons, ncol = 5))
colnames(panaroo_coincidence_binomial) <- c("gene1", "gene2", "obs_prop", "exp_prop", "p")

row_i <- 1

for (i in 1:length(panaroo_coincidence_binomial_raw)) {

  gene1 <- names(panaroo_coincidence_binomial_raw)[i]
  
  for(j in 1:length(panaroo_coincidence_binomial_raw[[i]])) {

    if (row_i %% 10000 == 0) { print(row_i) }
    
    gene2 <- names(panaroo_coincidence_binomial_raw[[i]])[j]

    panaroo_coincidence_binomial[row_i, c("gene1", "gene2")] <- c(gene1, gene2)

    panaroo_coincidence_binomial[row_i, c("obs_prop", "exp_prop", "p")] <- c(panaroo_coincidence_binomial_raw[[i]][[j]]$estimate,
                                                                             panaroo_coincidence_binomial_raw[[i]][[j]]$null.value,
                                                                             panaroo_coincidence_binomial_raw[[i]][[j]]$p.value)
    row_i <- row_i + 1

    }
}

saveRDS(object = panaroo_coincidence_binomial,
        file = "/home/gdouglas/projects/honey_bee/Ellegaard/coincidence_analysis/Gilliamella_and_Snodgrassella_panaroo_coincidence_binomial_df.rds")
