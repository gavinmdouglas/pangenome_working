rm(list = ls(all.names = TRUE))

library("discover")
library("dplyr")

Gilliamella_path_to_panaroo <- "/home/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Gilliamella_gene_presence_absence_roary-formatted.csv.gz"

Gilliamella_panaroo_out <- read.table(Gilliamella_path_to_panaroo,
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilli", rownames(Gilliamella_panaroo_out), sep = "_")

Gilliamella_panaroo_presence_raw <- Gilliamella_panaroo_out[, c(14:ncol(Gilliamella_panaroo_out))]
Gilliamella_panaroo_presence <- data.frame(matrix(NA, nrow = nrow(Gilliamella_panaroo_presence_raw), ncol = ncol(Gilliamella_panaroo_presence_raw)))
colnames(Gilliamella_panaroo_presence) <- colnames(Gilliamella_panaroo_presence_raw)
rownames(Gilliamella_panaroo_presence) <- rownames(Gilliamella_panaroo_presence_raw)
Gilliamella_panaroo_presence[Gilliamella_panaroo_presence_raw != ""] <- 1
Gilliamella_panaroo_presence[Gilliamella_panaroo_presence_raw == ""] <- 0

Gilliamella_ref_genome_events <- discover.matrix(Gilliamella_panaroo_presence)

Gilliamella_result.cooccur <- pairwise.discover.test(Gilliamella_ref_genome_events, fdr.method = "DBH", alternative = "greater")
Gilliamella_result.cooccur_df <- as.data.frame(Gilliamella_result.cooccur, q.threshold = 0.05)

Gilliamella_all_unique_hits <- c(Gilliamella_result.cooccur_df$gene1, Gilliamella_result.cooccur_df$gene2)
Gilliamella_all_unique_hits <- Gilliamella_all_unique_hits[-which(duplicated(Gilliamella_all_unique_hits))]

# Confirm that these genes co-occur in the actual data
Gilliamella_and_Gilliamella_present_noncore <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_and_Snodgrassella_present_noncore.rds")
Gilliamella_and_Gilliamella_present_noncore <- Gilliamella_and_Gilliamella_present_noncore[-which(rowSums(Gilliamella_and_Gilliamella_present_noncore) == 0), ]
Gilliamella_and_Gilliamella_present_noncore <- Gilliamella_and_Gilliamella_present_noncore[-which(rowSums(Gilliamella_and_Gilliamella_present_noncore) == 74), ]

presence_cutoff <- 0.5
Gilliamella_and_Gilliamella_present_noncore[Gilliamella_and_Gilliamella_present_noncore >= presence_cutoff] <- 1
Gilliamella_and_Gilliamella_present_noncore[Gilliamella_and_Gilliamella_present_noncore < presence_cutoff] <- 0

Ellegaard.2019_samples <- read.table("/home/gdouglas/projects/honey_bee/Ellegaard/SRR_ids.txt", stringsAsFactors = FALSE)$V1
Ellegaard.2020_samples <- read.table("/data1/gdouglas/projects/honey_bee/Ellegaard.2020/PRJNA598094_SRR_ids.txt", stringsAsFactors = FALSE)$V1

Gilliamella_and_Gilliamella_present_noncore <- Gilliamella_and_Gilliamella_present_noncore[, c(Ellegaard.2019_samples, Ellegaard.2020_samples)]

strata_vec <- c(rep("Switzerland", length(Ellegaard.2019_samples)), rep("Japan", length(Ellegaard.2020_samples)))

Gilliamella_breadth_events <- discover.matrix(Gilliamella_and_Gilliamella_present_noncore, strata = strata_vec)

Gilliamella_breadth_events_cooccur <- pairwise.discover.test(Gilliamella_breadth_events[rownames(Gilliamella_and_Gilliamella_present_noncore) %in% Gilliamella_all_unique_hits, ], fdr.method = "DBH", alternative = "greater")
Gilliamella_breadth_events_cooccur_df <- as.data.frame(Gilliamella_breadth_events_cooccur, q.threshold = 1)


Gilliamella_breadth_events_cooccur_df <- Gilliamella_breadth_events_cooccur_df %>%
  rowwise() %>%      # for each row
  mutate(compared = paste(sort(c(gene1, gene2)), collapse = " - ")) %>%  # sort the genes alphabetically and then combine them separating with -
  ungroup()          # forget the row grouping

Gilliamella_breadth_events_cooccur_df <- data.frame(Gilliamella_breadth_events_cooccur_df)

rownames(Gilliamella_breadth_events_cooccur_df) <- Gilliamella_breadth_events_cooccur_df$compared


Gilliamella_result.cooccur_df <- Gilliamella_result.cooccur_df %>%
  rowwise() %>%      # for each row
  mutate(compared = paste(sort(c(gene1, gene2)), collapse = " - ")) %>%  # sort the genes alphabetically and then combine them separating with -
  ungroup()          # forget the row grouping

Gilliamella_result.cooccur_df <- data.frame(Gilliamella_result.cooccur_df)

rownames(Gilliamella_result.cooccur_df) <- Gilliamella_result.cooccur_df$compared


# Only keep comparisons that were sig in the genome co-occurrence table.

Gilliamella_breadth_events_cooccur_df <- Gilliamella_breadth_events_cooccur_df[which(Gilliamella_breadth_events_cooccur_df$compared %in% Gilliamella_result.cooccur_df$compared), ]


Gilliamella_noncore_same_strain_cooccur <- list(genome_cooccur_results = Gilliamella_result.cooccur_df,
                                                  breadth_coccur_results = Gilliamella_breadth_events_cooccur_df)

saveRDS(object = Gilliamella_noncore_same_strain_cooccur,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Gilliamella_noncore_same_strain_cooccur.rds")



# Snodgrassella

rm(list = ls(all.names = TRUE))

Snodgrassella_path_to_panaroo <- "/home/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Snodgrassella_gene_presence_absence_roary-formatted.csv.gz"

Snodgrassella_panaroo_out <- read.table(Snodgrassella_path_to_panaroo,
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snod", rownames(Snodgrassella_panaroo_out), sep = "_")

Snodgrassella_panaroo_presence_raw <- Snodgrassella_panaroo_out[, c(14:ncol(Snodgrassella_panaroo_out))]
Snodgrassella_panaroo_presence <- data.frame(matrix(NA, nrow = nrow(Snodgrassella_panaroo_presence_raw), ncol = ncol(Snodgrassella_panaroo_presence_raw)))
colnames(Snodgrassella_panaroo_presence) <- colnames(Snodgrassella_panaroo_presence_raw)
rownames(Snodgrassella_panaroo_presence) <- rownames(Snodgrassella_panaroo_presence_raw)
Snodgrassella_panaroo_presence[Snodgrassella_panaroo_presence_raw != ""] <- 1
Snodgrassella_panaroo_presence[Snodgrassella_panaroo_presence_raw == ""] <- 0

Snodgrassella_ref_genome_events <- discover.matrix(Snodgrassella_panaroo_presence)

Snodgrassella_result.cooccur <- pairwise.discover.test(Snodgrassella_ref_genome_events, fdr.method = "DBH", alternative = "greater")
Snodgrassella_result.cooccur_df <- as.data.frame(Snodgrassella_result.cooccur, q.threshold = 0.05)

Snodgrassella_all_unique_hits <- c(Snodgrassella_result.cooccur_df$gene1, Snodgrassella_result.cooccur_df$gene2)
Snodgrassella_all_unique_hits <- Snodgrassella_all_unique_hits[-which(duplicated(Snodgrassella_all_unique_hits))]

# Confirm that these genes co-occur in the actual data
Snodgrassella_and_Snodgrassella_present_noncore <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_and_Snodgrassella_present_noncore.rds")
Snodgrassella_and_Snodgrassella_present_noncore <- Snodgrassella_and_Snodgrassella_present_noncore[-which(rowSums(Snodgrassella_and_Snodgrassella_present_noncore) == 0), ]
Snodgrassella_and_Snodgrassella_present_noncore <- Snodgrassella_and_Snodgrassella_present_noncore[-which(rowSums(Snodgrassella_and_Snodgrassella_present_noncore) == 74), ]

presence_cutoff <- 0.5
Snodgrassella_and_Snodgrassella_present_noncore[Snodgrassella_and_Snodgrassella_present_noncore >= presence_cutoff] <- 1
Snodgrassella_and_Snodgrassella_present_noncore[Snodgrassella_and_Snodgrassella_present_noncore < presence_cutoff] <- 0

Ellegaard.2019_samples <- read.table("/home/gdouglas/projects/honey_bee/Ellegaard/SRR_ids.txt", stringsAsFactors = FALSE)$V1
Ellegaard.2020_samples <- read.table("/data1/gdouglas/projects/honey_bee/Ellegaard.2020/PRJNA598094_SRR_ids.txt", stringsAsFactors = FALSE)$V1

Snodgrassella_and_Snodgrassella_present_noncore <- Snodgrassella_and_Snodgrassella_present_noncore[, c(Ellegaard.2019_samples, Ellegaard.2020_samples)]

strata_vec <- c(rep("Switzerland", length(Ellegaard.2019_samples)), rep("Japan", length(Ellegaard.2020_samples)))

Snodgrassella_breadth_events <- discover.matrix(Snodgrassella_and_Snodgrassella_present_noncore, strata = strata_vec)

Snodgrassella_breadth_events_cooccur <- pairwise.discover.test(Snodgrassella_breadth_events[rownames(Snodgrassella_and_Snodgrassella_present_noncore) %in% Snodgrassella_all_unique_hits, ], fdr.method = "DBH", alternative = "greater")
Snodgrassella_breadth_events_cooccur_df <- as.data.frame(Snodgrassella_breadth_events_cooccur, q.threshold = 1)


Snodgrassella_breadth_events_cooccur_df <- Snodgrassella_breadth_events_cooccur_df %>%
                                              rowwise() %>%      # for each row
                                              mutate(compared = paste(sort(c(gene1, gene2)), collapse = " - ")) %>%  # sort the genes alphabetically and then combine them separating with -
                                              ungroup()          # forget the row grouping

Snodgrassella_breadth_events_cooccur_df <- data.frame(Snodgrassella_breadth_events_cooccur_df)

rownames(Snodgrassella_breadth_events_cooccur_df) <- Snodgrassella_breadth_events_cooccur_df$compared


Snodgrassella_result.cooccur_df <- Snodgrassella_result.cooccur_df %>%
  rowwise() %>%      # for each row
  mutate(compared = paste(sort(c(gene1, gene2)), collapse = " - ")) %>%  # sort the genes alphabetically and then combine them separating with -
  ungroup()          # forget the row grouping

Snodgrassella_result.cooccur_df <- data.frame(Snodgrassella_result.cooccur_df)

rownames(Snodgrassella_result.cooccur_df) <- Snodgrassella_result.cooccur_df$compared


# Only keep comparisons that were sig in the genome co-occurrence table.
Snodgrassella_breadth_events_cooccur_df <- Snodgrassella_breadth_events_cooccur_df[which(Snodgrassella_breadth_events_cooccur_df$compared %in% Snodgrassella_result.cooccur_df$compared), ]


Snodgrassella_noncore_same_strain_cooccur <- list(genome_cooccur_results = Snodgrassella_result.cooccur_df,
                                                  breadth_coccur_results = Snodgrassella_breadth_events_cooccur_df)

saveRDS(object = Snodgrassella_noncore_same_strain_cooccur,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Snodgrassella_noncore_same_strain_cooccur.rds")

