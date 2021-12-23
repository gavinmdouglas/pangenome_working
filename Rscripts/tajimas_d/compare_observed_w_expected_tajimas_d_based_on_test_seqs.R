rm(list = ls(all.names = TRUE))

expected_popgen_df <- read.table(file = "honey_bee_pangenome/data/tests/test_seqs_for_tajimas_d/expected_pop_gen_metrics.tsv",
                                 header = TRUE, sep ="\t", row.names = 1)

multiplied_test <- read.table(file = "honey_bee_pangenome/data/tests/test_seqs_for_tajimas_d/sim_reads_popgen_pi_multiplied.tsv",
                              header = TRUE, sep ="\t", stringsAsFactors = FALSE, row.names = 1)

pi.per.site_test <- read.table(file = "honey_bee_pangenome/data/tests/test_seqs_for_tajimas_d/sim_reads_popgen_pi_per.site.tsv",
                              header = TRUE, sep ="\t", stringsAsFactors = FALSE, row.names = 1)

multiplied_instrain_test <- read.table(file = "honey_bee_pangenome/data/tests/test_seqs_for_tajimas_d/sim_reads_popgen_pi_instrain.tsv",
                              header = TRUE, sep ="\t", stringsAsFactors = FALSE, row.names = 1)

pi.per.site_instrain_test <- read.table(file = "honey_bee_pangenome/data/tests/test_seqs_for_tajimas_d/sim_reads_popgen_pi_instrain_per_site.tsv",
                                       header = TRUE, sep ="\t", stringsAsFactors = FALSE, row.names = 1)

rownames(multiplied_test) <- gsub("starting_seq", "rep", rownames(multiplied_test))
rownames(pi.per.site_test) <- gsub("starting_seq", "rep", rownames(pi.per.site_test))
rownames(multiplied_instrain_test) <- gsub("starting_seq", "rep", rownames(multiplied_instrain_test))
rownames(pi.per.site_instrain_test) <- gsub("starting_seq", "rep", rownames(pi.per.site_instrain_test))

multiplied_test <- multiplied_test[rownames(expected_popgen_df), ]
pi.per.site_test <- pi.per.site_test[rownames(expected_popgen_df), ]
multiplied_instrain_test <- multiplied_instrain_test[rownames(expected_popgen_df), ]
pi.per.site_instrain_test <- pi.per.site_instrain_test[rownames(expected_popgen_df), ]

plot(multiplied_test$tajimas_d, expected_popgen_df$Tajimas_D)
abline(a = 0, b = 1)

plot(multiplied_test$wattersons_theta, expected_popgen_df$Wattersons, xlim=c(0, 45), ylim=c(0, 45))
abline(a = 0, b = 1)

plot(multiplied_test$theta_pi, expected_popgen_df$abs_pi, xlim=c(0, 55), ylim=c(0, 55))
abline(a = 0, b = 1)

plot(multiplied_test$num_segregating_sites, expected_popgen_df$S)
abline(a = 0, b = 1)



plot(pi.per.site_test$wattersons_theta, expected_popgen_df$Wattersons)
abline(a = 0, b = 1)

plot(pi.per.site_test$theta_pi, expected_popgen_df$abs_pi, ylim=c(0, 55))
abline(a = 0, b = 1)

plot(pi.per.site_test$tajimas_d, expected_popgen_df$Tajimas_D)

plot(pi.per.site_test$num_segregating_sites, expected_popgen_df$S)


plot(multiplied_instrain_test$wattersons_theta, expected_popgen_df$Wattersons, xlim=c(0, 45), ylim=c(0, 45))
abline(a = 0, b = 1)

plot(multiplied_instrain_test$theta_pi, expected_popgen_df$abs_pi, ylim=c(0, 55))
abline(a = 0, b = 1)

plot(multiplied_instrain_test$tajimas_d, expected_popgen_df$Tajimas_D)
hist(multiplied_instrain_test$tajimas_d - expected_popgen_df$Tajimas_D)

plot(pi.per.site_instrain_test$wattersons_theta, expected_popgen_df$Wattersons, xlim=c(0, 45), ylim=c(0, 45))
abline(a = 0, b = 1)

plot(pi.per.site_instrain_test$theta_pi, expected_popgen_df$abs_pi, ylim=c(0, 55))
abline(a = 0, b = 1)

plot(pi.per.site_instrain_test$tajimas_d, expected_popgen_df$Tajimas_D)


plot(multiplied_instrain_test$num_segregating_sites - expected_popgen_df$S, multiplied_instrain_test$tajimas_d)

plot(multiplied_instrain_test$num_segregating_sites, multiplied_instrain_test$wattersons_theta)
abline(a = 0, b = 1)


library(ggplot2)

ggplot(data=multiplied_instrain_test, aes(x = tajimas_d, y = mean_polymorphic_coverage, colour = num_segregating_sites)) +
  geom_point()

plot(multiplied_instrain_test$tajimas_d, multiplied_instrain_test$mean_polymorphic_coverage)

calc_D_sd <- function(n, num_seg_sites) {
 
  a1 = sum(1 / c(1:(n - 1)))
  
  a2 = sum(1 / c(1:(n - 1))**2)
  
  b1 = (n + 1) / (3 * (n - 1))
  
  b2 = (2 * (n**2 + n + 3)) / (9 * n * (n - 1))
  
  c1 = b1 - (1 / a1)
  
  c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / a1**2)
  
  e1 = c1 / a1
  
  e2 = c2 / (a1**2 + a2)
  
  expected_sd = sqrt(e1 * num_seg_sites + e2 * num_seg_sites * (num_seg_sites - 1)) 
  
  return(expected_sd)
  
}


calc_a1 <- function(n) {
  
  return(sum(1 / c(1:(n - 1))))
  
}

multiplied_instrain_test$exp_std <- NA
expected_popgen_df$exp_std <- NA

multiplied_instrain_test$a1 <- NA
expected_popgen_df$a1 <- NA

for(i in 1:nrow(multiplied_instrain_test)) {
  multiplied_instrain_test[i, "exp_std"] <- calc_D_sd(n = multiplied_instrain_test[i, "mean_polymorphic_coverage"],
                                                      num_seg_sites = multiplied_instrain_test[i, "num_segregating_sites"])
  
  multiplied_instrain_test[i, "a1"] <- calc_a1(n = multiplied_instrain_test[i, "mean_polymorphic_coverage"])

  expected_popgen_df[i, "exp_std"] <- calc_D_sd(n = 10,
                                                num_seg_sites = expected_popgen_df[i, "S"])
  
  expected_popgen_df[i, "a1"] <- calc_a1(n = 10)
}



plot(multiplied_instrain_test$theta_pi / multiplied_instrain_test$num_segregating_sites, expected_popgen_df$Tajimas_D)
