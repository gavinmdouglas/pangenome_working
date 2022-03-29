
possible_strain_haplotype_combos_rep <- possible_strain_haplotype_combos_split[[which.max(combo_loglik_vec)]]

combo_loglik_vec[which.max(combo_loglik_vec)]

combo_likelihood <- 0

for (gene_sample in gene_and_strain_samples) {
  
  exp_prob <- compute_haplotype_exp_prob_based_on_strains(haplotype_dist = possible_strain_haplotype_combos_rep,
                                                          strain_abun = as.numeric(phylotype_strains[gene_sample, ]),
                                                          num_hap = num_haplotypes)
  
  sample_combo_likelihood <- NULL
  
  try(sample_combo_likelihood <- log(XNomial::xmulti(obs = sample_obs_abun[[gene_sample]], expr = exp_prob, detail = 0)$pProb), silent = TRUE)
  
  if (is.null(sample_combo_likelihood)) {
    sample_combo_likelihood <- log(XNomial::xmonte(obs = sample_obs_abun[[gene_sample]], expr = exp_prob, detail = 0)$observedProb)
  }
  
  combo_likelihood <- combo_likelihood + sample_combo_likelihood
  
}