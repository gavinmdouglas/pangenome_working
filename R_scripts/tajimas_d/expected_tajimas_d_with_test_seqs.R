rm(list = ls(all.names = TRUE))

library("ape")
library("pegas")

all_sim_fastas = list.files("honey_bee_pangenome/data/tests/test_seqs_for_tajimas_d/final_seqs/", pattern = "fasta")

expected_popgen_df <- data.frame(matrix(NA, nrow = length(all_sim_fastas), ncol = 6))
colnames(expected_popgen_df) <- c("rep", "pi_per_site", "abs_pi", "S", "Wattersons", "Tajimas_D")

expected_popgen_df$rep <- gsub(".fasta", "", all_sim_fastas)
expected_popgen_df$rep <- gsub("final_", "rep_", expected_popgen_df$rep)
rownames(expected_popgen_df) <- all_sim_fastas

for(sim_rep in all_sim_fastas) {
  
  dna_in <- read.dna(paste("honey_bee_pangenome/data/tests/test_seqs_for_tajimas_d/final_seqs/", sim_rep, sep = ""), format = "fasta")
  
  expected_popgen_df[sim_rep, c("pi_per_site", "abs_pi", "S", "Wattersons", "Tajimas_D")] <- c(nuc.div(dna_in),
                                                                                              nuc.div(dna_in) * 1000,
                                                                                              length(seg.sites(dna_in)),
                                                                                              theta.s(dna_in),
                                                                                              tajima.test(dna_in)$D)
  
}

write.table(x = expected_popgen_df, file = "honey_bee_pangenome/data/tests/test_seqs_for_tajimas_d/expected_pop_gen_metrics.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
