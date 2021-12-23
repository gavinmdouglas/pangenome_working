rm(list = ls(all.names = TRUE))

library(data.table)

panaroo_coincidence_binomial_raw <- readRDS("honey_bee_pangenome/data/Ellegaard_pangenome_coincidence/Gilliamella_and_Snodgrassella_panaroo_coincidence_binomial.rds")

# Remove the last element as it is just an empty list.
panaroo_coincidence_binomial_raw[length(panaroo_coincidence_binomial_raw)] <- NULL

p_val <- c()

row_i <- 1

for (i in 1:length(panaroo_coincidence_binomial_raw)) {
  
  gene1 <- names(panaroo_coincidence_binomial_raw)[i]
  
  for(j in 1:length(panaroo_coincidence_binomial_raw[[i]])) {
    
    if (row_i %% 10000 == 0) { print(row_i) }
    
    gene2 <- names(panaroo_coincidence_binomial_raw[[i]])[j]
    
    p_val <- c(p_val, panaroo_coincidence_binomial_raw[[i]][[j]]$p.value)
    
    #    panaroo_coincidence_binomial <- rbind(panaroo_coincidence_binomial, list(gene1,
    #                                                                            gene2,
    #                                                                            panaroo_coincidence_binomial_raw[[i]][[j]]$estimate,
    #                                                                            panaroo_coincidence_binomial_raw[[i]][[j]]$null.value,
    #                                                                            panaroo_coincidence_binomial_raw[[i]][[j]]$p.value))
    
    row_i <- row_i + 1
    
  }
}
