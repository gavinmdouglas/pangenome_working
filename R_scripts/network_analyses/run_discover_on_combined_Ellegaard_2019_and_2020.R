### Code to run DISCOVER analysis between all genes  ###
### Run co-exclusion and co-occurence separately and also run 2019 and 2020 data together (stratified) ### 

rm(list = ls(all.names = TRUE))

library("discover")

pandora_output_noncore <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.noncore.matrix",
                             header = TRUE, sep = "\t", row.names = 1)

# Remove non-core genes present in all or missing in all samples.
pandora_output_noncore <- pandora_output_noncore[-which(rowSums(pandora_output_noncore) == ncol(pandora_output_noncore)), ]
#pandora_output_noncore <- pandora_output_noncore[-which(rowSums(pandora_output_noncore) == 0), ]

Ellegaard.2019_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2019_SRR_ids.txt", stringsAsFactors = FALSE)$V1
Ellegaard.2020_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2020_SRR_ids.txt", stringsAsFactors = FALSE)$V1

Ellegaard.2019_samples <- Ellegaard.2019_samples[which(Ellegaard.2019_samples %in% colnames(pandora_output_noncore))]
Ellegaard.2020_samples <- Ellegaard.2020_samples[which(Ellegaard.2020_samples %in% colnames(pandora_output_noncore))]

pandora_output_noncore <- pandora_output_noncore[, c(Ellegaard.2019_samples, Ellegaard.2020_samples)]

# Remove Bombilactobacillus_mellifer, Gilliamella_apis, Frischella_perrara, and Commensalibacter_sp
# as < 20 genes for each of those species were in this subset.
pandora_output_noncore <- pandora_output_noncore[-grep("Bombilactobacillus_mellifer", rownames(pandora_output_noncore)), ]
pandora_output_noncore <- pandora_output_noncore[-grep("Gilliamella_apis", rownames(pandora_output_noncore)), ]
pandora_output_noncore <- pandora_output_noncore[-grep("Frischella_perrara", rownames(pandora_output_noncore)), ]
pandora_output_noncore <- pandora_output_noncore[-grep("Commensalibacter_sp", rownames(pandora_output_noncore)), ]

strata_vec <- c(rep("Switzerland", length(Ellegaard.2019_samples)), rep("Japan", length(Ellegaard.2020_samples)))

events <- discover.matrix(pandora_output_noncore, strata = strata_vec)

subset_to_test <- rowSums(pandora_output_noncore) >= 10 & rowSums(pandora_output_noncore == 0) >= 10

saveRDS(object = subset_to_test,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_tested_genes.rds")

result.mutex <- pairwise.discover.test(events[subset_to_test, ], fdr.method = "DBH", alternative = "less")
result.mutex_df <- as.data.frame(result.mutex, q.threshold = 1)

result.cooccur <- pairwise.discover.test(events[subset_to_test, ], fdr.method = "DBH", alternative = "greater")
result.cooccur_df <- as.data.frame(result.cooccur, q.threshold = 1)

result.mutex_df_DBH0.35 <- result.mutex_df[which(result.mutex_df$q.value < 0.35), ]

result.cooccur_df_DBH0.35 <- result.cooccur_df[which(result.cooccur_df$q.value < 0.35), ]


saveRDS(object = result.mutex_df_DBH0.35,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_mutex_results_clustered_DBH0.35.rds")

saveRDS(object = result.cooccur_df_DBH0.35,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_cooccur_results_clustered_DBH0.35.rds")
