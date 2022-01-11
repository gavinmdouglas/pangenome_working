rm(list = ls(all.names = TRUE))


par(mfrow=c(1, 2))
Gilliamella_SRR7287194_test <- read.table("honey_bee_pangenome/data/Ellegaard_Carrie_pangenome_mapping/popgen/subfam80_test/Gilliamella/Gilliamella.SRR7287194.subfam80.popgen.tsv.gz",
                                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)

plot(Gilliamella_SRR7287194_test$theta_pi, Gilliamella_SRR7287194_test$wattersons_theta,
     main = "Gilliamella - SRR7287194", xlab = "Theta pi", ylab = "Watterson's theta")

hist(Gilliamella_SRR7287194_test$tajimas_d, main = "Gilliamella - SRR7287194", xlab = "Tajima's D")




Gilliamella_SRR7287230_test <- read.table("honey_bee_pangenome/data/Ellegaard_Carrie_pangenome_mapping/popgen/subfam80_test/Gilliamella/Gilliamella.SRR7287230.subfam80.popgen.tsv.gz",
                                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
plot(Gilliamella_SRR7287230_test$theta_pi, Gilliamella_SRR7287194_test$wattersons_theta,
     main = "Gilliamella - SRR7287230", xlab = "Theta pi", ylab = "Watterson's theta")

hist(Gilliamella_SRR7287230_test$tajimas_d, main = "Gilliamella - SRR7287230", xlab = "Tajima's D")



Gilliamella_SRR7287194_test <- read.table("honey_bee_pangenome/data/Ellegaard_Carrie_pangenome_mapping/popgen/subfam80_test/Gilliamella/Gilliamella.SRR7287194.subfam80.popgen.tsv.gz",
                                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
plot(Gilliamella_SRR7287194_test$theta_pi, Gilliamella_SRR7287194_test$wattersons_theta,
     main = "Gilliamella - SRR7287194", xlab = "Theta pi", ylab = "Watterson's theta")

hist(Gilliamella_SRR7287194_test$tajimas_d, main = "Gilliamella - SRR7287194", xlab = "Tajima's D")


tajimas_d_sum <- rep(0, nrow(Gilliamella_SRR7287194_test))
tajimas_nonNA_sum <- rep(0, nrow(Gilliamella_SRR7287194_test))

Gilliamella_pop_gen_files <- list.files("honey_bee_pangenome/data/Ellegaard_Carrie_pangenome_mapping/popgen/subfam80_test/Gilliamella/",
                                        full.names = TRUE)

for(f in Gilliamella_pop_gen_files) {
  
  tmp <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  tajimas_nonNA_sum <- tajimas_nonNA_sum + ! is.na(tmp$tajimas_d)
  
  tmp$tajimas_d[which(is.na(tmp$tajimas_d))] <- 0
  
  tajimas_d_sum <- tajimas_d_sum + tmp$tajimas_d
   
}

tajimas_d_mean <- tajimas_d_sum / tajimas_nonNA_sum
hist(tajimas_d_mean, main = "Gilliamella", xlab = "Mean Tajima's D (per gene)")
plot(tajimas_nonNA_sum, tajimas_d_mean, xlab = "No. samples", ylab = "Tajima's D", main = "Gilliamella")


tajimas_d_mean_char <- as.character(tajimas_d_mean)

0.2792417

tajimas_nonNA_sum_tm


par(mfrow=c(1, 2))
Snodgrassella_SRR7287194_test <- read.table("honey_bee_pangenome/data/Ellegaard_Carrie_pangenome_mapping/popgen/subfam80_test/Snodgrassella/Snodgrassella.SRR7287194.subfam80.popgen.tsv.gz",
                                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)

plot(Snodgrassella_SRR7287194_test$theta_pi, Snodgrassella_SRR7287194_test$wattersons_theta,
     main = "Snodgrassella - SRR7287194", xlab = "Theta pi", ylab = "Watterson's theta")

hist(Snodgrassella_SRR7287194_test$tajimas_d, main = "Snodgrassella - SRR7287194", xlab = "Tajima's D")




Snodgrassella_SRR7287230_test <- read.table("honey_bee_pangenome/data/Ellegaard_Carrie_pangenome_mapping/popgen/subfam80_test/Snodgrassella/Snodgrassella.SRR7287230.subfam80.popgen.tsv.gz",
                                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
plot(Snodgrassella_SRR7287230_test$theta_pi, Snodgrassella_SRR7287194_test$wattersons_theta,
     main = "Snodgrassella - SRR7287230", xlab = "Theta pi", ylab = "Watterson's theta")

hist(Snodgrassella_SRR7287230_test$tajimas_d, main = "Snodgrassella - SRR7287230", xlab = "Tajima's D")



Snodgrassella_SRR7287194_test <- read.table("honey_bee_pangenome/data/Ellegaard_Carrie_pangenome_mapping/popgen/subfam80_test/Snodgrassella/Snodgrassella.SRR7287194.subfam80.popgen.tsv.gz",
                                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)

plot(Snodgrassella_SRR7287194_test$theta_pi, Snodgrassella_SRR7287194_test$wattersons_theta,
       main = "Snodgrassella - SRR7287194", xlab = "Theta pi", ylab = "Watterson's theta")

hist(Snodgrassella_SRR7287194_test$tajimas_d, main = "Snodgrassella - SRR7287194", xlab = "Tajima's D")


tajimas_d_sum <- rep(0, nrow(Snodgrassella_SRR7287194_test))
tajimas_nonNA_sum <- rep(0, nrow(Snodgrassella_SRR7287194_test))

Snodgrassella_pop_gen_files <- list.files("honey_bee_pangenome/data/Ellegaard_Carrie_pangenome_mapping/popgen/subfam80_test/Snodgrassella/",
                                        full.names = TRUE)

for(f in Snodgrassella_pop_gen_files) {
  
  tmp <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  tajimas_nonNA_sum <- tajimas_nonNA_sum + ! is.na(tmp$tajimas_d)
  
  tmp$tajimas_d[which(is.na(tmp$tajimas_d))] <- 0
  
  tajimas_d_sum <- tajimas_d_sum + tmp$tajimas_d
  
}

tajimas_d_mean <- tajimas_d_sum / tajimas_nonNA_sum
hist(tajimas_d_mean, main = "Snodgrassella", xlab = "Mean Tajima's D (per gene)")
plot(tajimas_nonNA_sum, tajimas_d_mean, xlab = "No. samples", ylab = "Tajima's D", main = "Snodgrassella")