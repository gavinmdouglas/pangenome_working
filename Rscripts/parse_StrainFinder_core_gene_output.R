# Parse AIC files for all StrainFinder core gene outputs and then make a list with the OTU table for each phylotype.
# Need to add in the sample ids to this table.
# Save this as an RDS so that the data can be quickly reused in the future.

rm(list = ls(all.names = TRUE))

# Chose this cut-off based on distribution of strains across samples (there tended to just be major outlier after this point in most cases)
# Also, this corresponds to a relative abundance of 0.01% - it doesn't seem reasonable to expect to be able to identify strains at lower abundance than that.
min_abun <- 0.01

abun_tables <- list()

Bifidobacterium_pandora.struct_AICs <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_genes_output/pandora_prepped_w_struct/Bifidobacterium_core.strain_fit_summary.tsv",
                                                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)


Bifidobacterium_pandora.struct_num_strain_best_fit <- as.character(Bifidobacterium_pandora.struct_AICs[which(Bifidobacterium_pandora.struct_AICs$AIC == min(Bifidobacterium_pandora.struct_AICs$AIC)), "ID"])

Bifidobacterium_pandora.struct_num_strain_best_fit_abun_path <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_genes_output/pandora_prepped_w_struct/otu_table.Bifidobacterium_core.",
                                                                      Bifidobacterium_pandora.struct_num_strain_best_fit,
                                                                      ".txt",
                                                                      sep = "")

Bifidobacterium_pandora.struct_num_strain_best_fit_abun <- read.table(file = Bifidobacterium_pandora.struct_num_strain_best_fit_abun_path,
                                                                      header = FALSE, sep = "\t", skip = 1)

Bifidobacterium_core_pandora.struct_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/pandora_prepped_core_gene_StrainFinder_input/prepped_input/Bifidobacterium_core_samples.txt",
                                                          header = FALSE, sep = "", stringsAsFactors = FALSE)$V1

rownames(Bifidobacterium_pandora.struct_num_strain_best_fit_abun) <- Bifidobacterium_core_pandora.struct_samples

Bifidobacterium_pandora.struct_num_strain_best_fit_abun_filt <- Bifidobacterium_pandora.struct_num_strain_best_fit_abun
Bifidobacterium_pandora.struct_num_strain_best_fit_abun_filt[Bifidobacterium_pandora.struct_num_strain_best_fit_abun_filt < min_abun] <- 0
Bifidobacterium_pandora.struct_num_strain_best_fit_abun_filt_relabun <- data.frame(sweep(x = Bifidobacterium_pandora.struct_num_strain_best_fit_abun_filt,
                                                                                         MARGIN = 1,
                                                                                         STATS = rowSums(Bifidobacterium_pandora.struct_num_strain_best_fit_abun_filt),
                                                                                         FUN = '/')) * 100

abun_tables[["Bifidobacterium"]] <- Bifidobacterium_pandora.struct_num_strain_best_fit_abun_filt_relabun






Firm4_pandora.struct_AICs <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_genes_output/pandora_prepped_w_struct//Firm4_core.strain_fit_summary.tsv",
                                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)


Firm4_pandora.struct_num_strain_best_fit <- as.character(Firm4_pandora.struct_AICs[which(Firm4_pandora.struct_AICs$AIC == min(Firm4_pandora.struct_AICs$AIC)), "ID"])

Firm4_pandora.struct_num_strain_best_fit_abun_path <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_genes_output/pandora_prepped_w_struct/otu_table.Firm4_core.",
                                                                      Firm4_pandora.struct_num_strain_best_fit,
                                                                      ".txt",
                                                                      sep = "")

Firm4_pandora.struct_num_strain_best_fit_abun <- read.table(file = Firm4_pandora.struct_num_strain_best_fit_abun_path,
                                                                      header = FALSE, sep = "\t", skip = 1)

Firm4_core_pandora.struct_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/pandora_prepped_core_gene_StrainFinder_input/prepped_input/Firm4_core_samples.txt",
                                                          header = FALSE, sep = "", stringsAsFactors = FALSE)$V1

rownames(Firm4_pandora.struct_num_strain_best_fit_abun) <- Firm4_core_pandora.struct_samples

Firm4_pandora.struct_num_strain_best_fit_abun_filt <- Firm4_pandora.struct_num_strain_best_fit_abun
Firm4_pandora.struct_num_strain_best_fit_abun_filt[Firm4_pandora.struct_num_strain_best_fit_abun_filt < min_abun] <- 0
Firm4_pandora.struct_num_strain_best_fit_abun_filt_relabun <- data.frame(sweep(x = Firm4_pandora.struct_num_strain_best_fit_abun_filt,
                                                                                         MARGIN = 1,
                                                                                         STATS = rowSums(Firm4_pandora.struct_num_strain_best_fit_abun_filt),
                                                                                         FUN = '/')) * 100

abun_tables[["Firm4"]] <- Firm4_pandora.struct_num_strain_best_fit_abun_filt_relabun









Firm5_pandora.struct_AICs <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_genes_output/pandora_prepped_w_struct//Firm5_core.strain_fit_summary.tsv",
                                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)


Firm5_pandora.struct_num_strain_best_fit <- as.character(Firm5_pandora.struct_AICs[which(Firm5_pandora.struct_AICs$AIC == min(Firm5_pandora.struct_AICs$AIC)), "ID"])

Firm5_pandora.struct_num_strain_best_fit_abun_path <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_genes_output/pandora_prepped_w_struct/otu_table.Firm5_core.",
                                                            Firm5_pandora.struct_num_strain_best_fit,
                                                            ".txt",
                                                            sep = "")

Firm5_pandora.struct_num_strain_best_fit_abun <- read.table(file = Firm5_pandora.struct_num_strain_best_fit_abun_path,
                                                            header = FALSE, sep = "\t", skip = 1)

Firm5_core_pandora.struct_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/pandora_prepped_core_gene_StrainFinder_input/prepped_input/Firm5_core_samples.txt",
                                                header = FALSE, sep = "", stringsAsFactors = FALSE)$V1

rownames(Firm5_pandora.struct_num_strain_best_fit_abun) <- Firm5_core_pandora.struct_samples

Firm5_pandora.struct_num_strain_best_fit_abun_filt <- Firm5_pandora.struct_num_strain_best_fit_abun
Firm5_pandora.struct_num_strain_best_fit_abun_filt[Firm5_pandora.struct_num_strain_best_fit_abun_filt < min_abun] <- 0
Firm5_pandora.struct_num_strain_best_fit_abun_filt_relabun <- data.frame(sweep(x = Firm5_pandora.struct_num_strain_best_fit_abun_filt,
                                                                               MARGIN = 1,
                                                                               STATS = rowSums(Firm5_pandora.struct_num_strain_best_fit_abun_filt),
                                                                               FUN = '/')) * 100

abun_tables[["Firm5"]] <- Firm5_pandora.struct_num_strain_best_fit_abun_filt_relabun







Gilliamella_pandora.struct_AICs <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_genes_output/pandora_prepped_w_struct//Gilliamella_core.strain_fit_summary.tsv",
                                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)


Gilliamella_pandora.struct_num_strain_best_fit <- as.character(Gilliamella_pandora.struct_AICs[which(Gilliamella_pandora.struct_AICs$AIC == min(Gilliamella_pandora.struct_AICs$AIC)), "ID"])

Gilliamella_pandora.struct_num_strain_best_fit_abun_path <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_genes_output/pandora_prepped_w_struct/otu_table.Gilliamella_core.",
                                                            Gilliamella_pandora.struct_num_strain_best_fit,
                                                            ".txt",
                                                            sep = "")

Gilliamella_pandora.struct_num_strain_best_fit_abun <- read.table(file = Gilliamella_pandora.struct_num_strain_best_fit_abun_path,
                                                            header = FALSE, sep = "\t", skip = 1)

Gilliamella_core_pandora.struct_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/pandora_prepped_core_gene_StrainFinder_input/prepped_input/Gilliamella_core_samples.txt",
                                                header = FALSE, sep = "", stringsAsFactors = FALSE)$V1

rownames(Gilliamella_pandora.struct_num_strain_best_fit_abun) <- Gilliamella_core_pandora.struct_samples

Gilliamella_pandora.struct_num_strain_best_fit_abun_filt <- Gilliamella_pandora.struct_num_strain_best_fit_abun
Gilliamella_pandora.struct_num_strain_best_fit_abun_filt[Gilliamella_pandora.struct_num_strain_best_fit_abun_filt < min_abun] <- 0
Gilliamella_pandora.struct_num_strain_best_fit_abun_filt_relabun <- data.frame(sweep(x = Gilliamella_pandora.struct_num_strain_best_fit_abun_filt,
                                                                               MARGIN = 1,
                                                                               STATS = rowSums(Gilliamella_pandora.struct_num_strain_best_fit_abun_filt),
                                                                               FUN = '/')) * 100

abun_tables[["Gilliamella"]] <- Gilliamella_pandora.struct_num_strain_best_fit_abun_filt_relabun








Snodgrassella_pandora.struct_AICs <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_genes_output/pandora_prepped_w_struct//Snodgrassella_core.strain_fit_summary.tsv",
                                                header = TRUE, sep = "\t", stringsAsFactors = FALSE)


Snodgrassella_pandora.struct_num_strain_best_fit <- as.character(Snodgrassella_pandora.struct_AICs[which(Snodgrassella_pandora.struct_AICs$AIC == min(Snodgrassella_pandora.struct_AICs$AIC)), "ID"])

Snodgrassella_pandora.struct_num_strain_best_fit_abun_path <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_genes_output/pandora_prepped_w_struct/otu_table.Snodgrassella_core.",
                                                                    Snodgrassella_pandora.struct_num_strain_best_fit,
                                                                    ".txt",
                                                                    sep = "")

Snodgrassella_pandora.struct_num_strain_best_fit_abun <- read.table(file = Snodgrassella_pandora.struct_num_strain_best_fit_abun_path,
                                                                    header = FALSE, sep = "\t", skip = 1)

Snodgrassella_core_pandora.struct_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/pandora_prepped_core_gene_StrainFinder_input/prepped_input/Snodgrassella_core_samples.txt",
                                                        header = FALSE, sep = "", stringsAsFactors = FALSE)$V1

rownames(Snodgrassella_pandora.struct_num_strain_best_fit_abun) <- Snodgrassella_core_pandora.struct_samples

Snodgrassella_pandora.struct_num_strain_best_fit_abun_filt <- Snodgrassella_pandora.struct_num_strain_best_fit_abun
Snodgrassella_pandora.struct_num_strain_best_fit_abun_filt[Snodgrassella_pandora.struct_num_strain_best_fit_abun_filt < min_abun] <- 0
Snodgrassella_pandora.struct_num_strain_best_fit_abun_filt_relabun <- data.frame(sweep(x = Snodgrassella_pandora.struct_num_strain_best_fit_abun_filt,
                                                                                     MARGIN = 1,
                                                                                     STATS = rowSums(Snodgrassella_pandora.struct_num_strain_best_fit_abun_filt),
                                                                                     FUN = '/')) * 100

abun_tables[["Snodgrassella"]] <- Snodgrassella_pandora.struct_num_strain_best_fit_abun_filt_relabun

saveRDS(object = abun_tables,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/core_genes_output/strainfinder_pandora_struct_abun.rds")
