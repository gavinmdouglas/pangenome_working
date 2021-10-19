rm(list = ls(all.names = TRUE))

library(data.table)
library(parallel)

panaroo_coincidence_binomial_raw <- readRDS("/home/gdouglas/projects/honey_bee/Ellegaard/coincidence_analysis/Gilliamella_and_Snodgrassella_panaroo_coincidence_binomial.rds")
#panaroo_coincidence_binomial_raw <- readRDS("honey_bee_pangenome/data/Ellegaard_pangenome_coincidence/Gilliamella_and_Snodgrassella_panaroo_coincidence_binomial.rds")

# Remove the last element as it is just an empty list.
panaroo_coincidence_binomial_raw[length(panaroo_coincidence_binomial_raw)] <- NULL

parse_binomial_test_info <- function(binom_test_output) {
 
  out_df <- data.table(gene1 = character(),
                       gene2 = character(),
                       obs_prop = numeric(),
                       exp_prop = numeric(),
                       p = numeric())
   
  for(j in 1:length(binom_test_output)) {
    
    gene2 <- names(binom_test_output)[j]
  
    out_df <- rbind(out_df, list(NA,
            								     gene2,
            								     binom_test_output[[j]]$estimate,
            								     binom_test_output[[j]]$null.value,
            								     binom_test_output[[j]]$p.value))
    
  }
  
  return(out_df)
  
}

panaroo_coincidence_binomial_raw_parsed <- mclapply(X = panaroo_coincidence_binomial_raw, FUN = parse_binomial_test_info, mc.cores = 30)


for (i in 1:length(panaroo_coincidence_binomial_raw_parsed)) {
  panaroo_coincidence_binomial_raw_parsed[[i]]$gene1 <- names(panaroo_coincidence_binomial_raw_parsed)[i]
}

panaroo_coincidence_binomial <- do.call(rbind, panaroo_coincidence_binomial_raw_parsed)

panaroo_coincidence_binomial$BH <- p.adjust(panaroo_coincidence_binomial$p, "BH")

saveRDS(object = panaroo_coincidence_binomial,
        file = "/home/gdouglas/projects/honey_bee/Ellegaard/coincidence_analysis/Gilliamella_and_Snodgrassella_panaroo_coincidence_binomial_df.rds")

panaroo_coincidence_binomial_BH_0.4 <- panaroo_coincidence_binomial[which(panaroo_coincidence_binomial$BH < 0.4), ]

saveRDS(object = panaroo_coincidence_binomial_BH_0.4,
        file = "/home/gdouglas/projects/honey_bee/Ellegaard/coincidence_analysis/Gilliamella_and_Snodgrassella_panaroo_coincidence_binomial_BH_0.4_df.rds")



panaroo_coincidence_binomial_Gilli <- data.frame(panaroo_coincidence_binomial)
panaroo_coincidence_binomial_Gilli <- panaroo_coincidence_binomial_Gilli[-grep("Snod", panaroo_coincidence_binomial_Gilli$gene1), ]
panaroo_coincidence_binomial_Gilli <- panaroo_coincidence_binomial_Gilli[-grep("Snod", panaroo_coincidence_binomial_Gilli$gene2), ]


panaroo_coincidence_binomial_Snod <- data.frame(panaroo_coincidence_binomial)
panaroo_coincidence_binomial_Snod <- panaroo_coincidence_binomial_Snod[-grep("Gilli", panaroo_coincidence_binomial_Snod$gene1), ]
