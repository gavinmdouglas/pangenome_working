### Code to run DISCOVER analysis between all Gilliamella and Snodgrassella genes within a certain window. ###
### Run coexclusion and cooccurence separately and also run 2019 and 2020 data together (stratified) ### 

rm(list = ls(all.names = TRUE))

library("discover")

Gilliamella_and_Snodgrassella_present_noncore <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_and_Snodgrassella_present_noncore.rds")

Gilliamella_and_Snodgrassella_present_noncore <- Gilliamella_and_Snodgrassella_present_noncore[-which(rowSums(Gilliamella_and_Snodgrassella_present_noncore) == 0), ]

presence_cutoff <- 0.5
Gilliamella_and_Snodgrassella_present_noncore[Gilliamella_and_Snodgrassella_present_noncore >= presence_cutoff] <- 1
Gilliamella_and_Snodgrassella_present_noncore[Gilliamella_and_Snodgrassella_present_noncore < presence_cutoff] <- 0

Ellegaard.2019_samples <- read.table("/home/gdouglas/projects/honey_bee/Ellegaard/SRR_ids.txt", stringsAsFactors = FALSE)$V1
Ellegaard.2020_samples <- read.table("/data1/gdouglas/projects/honey_bee/Ellegaard.2020/PRJNA598094_SRR_ids.txt", stringsAsFactors = FALSE)$V1

Gilliamella_and_Snodgrassella_present_noncore <- Gilliamella_and_Snodgrassella_present_noncore[, c(Ellegaard.2019_samples, Ellegaard.2020_samples)]

strata_vec <- c(rep("Switzerland", length(Ellegaard.2019_samples)), rep("Japan", length(Ellegaard.2020_samples)))

events <- discover.matrix(Gilliamella_and_Snodgrassella_present_noncore, strata = strata_vec)

# Note that it's important that to subset after running discover.matrix
subset_to_test <- rowSums(Gilliamella_and_Snodgrassella_present_noncore) >= 5 & rowSums(Gilliamella_and_Snodgrassella_present_noncore) <= 0.75 * ncol(Gilliamella_and_Snodgrassella_present_noncore)

result.mutex <- pairwise.discover.test(events[subset_to_test, ], fdr.method = "DBH", alternative = "less")
result.mutex_df <- as.data.frame(result.mutex, q.threshold = 0.3)

result.cooccur <- pairwise.discover.test(events[subset_to_test, ], fdr.method = "DBH", alternative = "greater")
result.cooccur_df <- as.data.frame(result.cooccur, q.threshold = 0.3)

Gilli_gene1 <- grep("Gilli", result.cooccur_df$gene1)
Gilli_gene2 <- grep("Gilli", result.cooccur_df$gene2)
Snod_gene1 <- grep("Snod", result.cooccur_df$gene1)
Snod_gene2 <- grep("Snod", result.cooccur_df$gene2)

length(which(Gilli_gene1 %in% Snod_gene2)) + length(which(Gilli_gene2 %in% Snod_gene1))

saveRDS(object = result.cooccur_df,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/DISCOVER_network/Ellegaard_2019_2020_DISCOVER_cooccur_results_DBH0.3.rds")
