### Code to run DISCOVER analysis between all Gilliamella and Snodgrassella genes within a certain window. ###
### Run co-exclusion and co-occurence separately and also run 2019 and 2020 data together (stratified) ### 

rm(list = ls(all.names = TRUE))

library("discover")

pandora_output_noncore <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_multisample.filt.noncore.matrix",
                             header = TRUE, sep = "\t", row.names = 1)

# Remove non-core genes present in all or missing in all samples.
pandora_output_noncore <- pandora_output_noncore[-which(rowSums(pandora_output_noncore) == ncol(pandora_output_noncore)), ]
pandora_output_noncore <- pandora_output_noncore[-which(rowSums(pandora_output_noncore) == 0), ]

Ellegaard.2019_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2019_SRR_ids.txt", stringsAsFactors = FALSE)$V1
Ellegaard.2020_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2020_SRR_ids.txt", stringsAsFactors = FALSE)$V1

Ellegaard.2019_samples <- Ellegaard.2019_samples[which(Ellegaard.2019_samples %in% colnames(pandora_output_noncore))]
Ellegaard.2020_samples <- Ellegaard.2020_samples[which(Ellegaard.2020_samples %in% colnames(pandora_output_noncore))]

pandora_output_noncore <- pandora_output_noncore[, c(Ellegaard.2019_samples, Ellegaard.2020_samples)]

strata_vec <- c(rep("Switzerland", length(Ellegaard.2019_samples)), rep("Japan", length(Ellegaard.2020_samples)))

events <- discover.matrix(pandora_output_noncore, strata = strata_vec)

subset_to_test <- rowSums(pandora_output_noncore) >= 5 & rowSums(pandora_output_noncore) <= 67

saveRDS(object = subset_to_test,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_tested_genes.rds")

result.mutex <- pairwise.discover.test(events[subset_to_test, ], fdr.method = "DBH", alternative = "less")
result.mutex_df <- as.data.frame(result.mutex, q.threshold = 1)

result.cooccur <- pairwise.discover.test(events[subset_to_test, ], fdr.method = "DBH", alternative = "greater")
result.cooccur_df <- as.data.frame(result.cooccur, q.threshold = 1)

result.cooccur_df_DBH0.3 <- result.cooccur_df[which(result.cooccur_df$q.value < 0.3), ]

saveRDS(object = result.cooccur_df_DBH0.3,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_cooccur_results_clustered_DBH0.3.rds")

