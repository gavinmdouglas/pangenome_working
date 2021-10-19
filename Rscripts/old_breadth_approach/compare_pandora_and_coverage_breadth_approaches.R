### Compare inferred presence / absence profile based on pandora vs coverage breadth approach

rm(list = ls(all.names = TRUE))

pandora_output <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_multisample.matrix",
                             header = TRUE, sep = "\t", row.names = 1)
rownames(pandora_output) <- gsub("\\.fa$", "", rownames(pandora_output))

breadth_out <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_and_Snodgrassella_present.rds")

breadth_out_0.1 <- breadth_out
breadth_out_0.5 <- breadth_out
breadth_out_0.9 <- breadth_out

breadth_out_0.1[breadth_out_0.1 < 0.1] <- 0
breadth_out_0.1[breadth_out_0.1 >= 0.1] <- 1
breadth_out_0.1 <- breadth_out_0.1[-which(rowSums(breadth_out_0.1) == 0), ]

breadth_out_0.5[breadth_out_0.5 < 0.5] <- 0
breadth_out_0.5[breadth_out_0.5 >= 0.5] <- 1
breadth_out_0.5 <- breadth_out_0.5[-which(rowSums(breadth_out_0.5) == 0), ]

breadth_out_0.9[breadth_out_0.9 < 0.9] <- 0
breadth_out_0.9[breadth_out_0.9 >= 0.9] <- 1
breadth_out_0.9 <- breadth_out_0.9[-which(rowSums(breadth_out_0.9) == 0), ]

length(grep("Gilli_", rownames(pandora_output)))
length(grep("Gilli_", rownames(breadth_out_0.1)))
length(grep("Gilli_", rownames(breadth_out_0.5)))
length(grep("Gilli_", rownames(breadth_out_0.9)))

intersecting_genes <- rownames(pandora_output)[which(rownames(pandora_output) %in% rownames(breadth_out_0.5))]
intersecting_genes_Gilli <- intersecting_genes[grep("Gilli_", intersecting_genes)]
intersecting_genes_Snod <- intersecting_genes[grep("Snod_", intersecting_genes)]

### Compare observed and expected intersect of positive calls

Gilli_obs <- c()
Gilli_exp <- c()

Snod_obs <- c()
Snod_exp <- c()

for (SRR in colnames(pandora_output)) {
  pandora_Gilli_present_calls <- which(pandora_output[intersecting_genes_Gilli, SRR] == 1)
  breadth_Gilli_present_calls <- which(breadth_out_0.5[intersecting_genes_Gilli, SRR] == 1)
  Gilli_obs <- c(Gilli_obs, length(which(pandora_Gilli_present_calls %in% breadth_Gilli_present_calls)))
  Gilli_exp <- c(Gilli_exp, length(pandora_Gilli_present_calls) * (length(breadth_Gilli_present_calls) / length(intersecting_genes_Gilli)))
  
  pandora_Snod_present_calls <- which(pandora_output[intersecting_genes_Snod, SRR] == 1)
  breadth_Snod_present_calls <- which(breadth_out_0.5[intersecting_genes_Snod, SRR] == 1)
  Snod_obs <- c(Snod_obs, length(which(pandora_Snod_present_calls %in% breadth_Snod_present_calls)))
  Snod_exp <- c(Snod_exp, length(pandora_Snod_present_calls) * (length(breadth_Snod_present_calls) / length(intersecting_genes_Snod)))
  
}

par(mfrow=c(2,1))
plot(Gilli_obs, Gilli_exp, main = "Gilliamella", xlab = "Observed intersecting positive calls", ylab = "Expected intersecting positive calls at random")
abline(b = 1, a = 0, lty = 2, lwd = 2)

plot(Snod_obs, Snod_exp, main = "Snodgrassella", xlab = "Observed intersecting positive calls", ylab = "Expected intersecting positive calls at random")
abline(b = 1, a = 0, lty = 2, lwd = 2)
par(mfrow=c(1,1))



### Compare total numbers of positive calls per sample.
Gilli_pandora_counts <- c()
Gilli_breadth_counts <- c()

Snod_pandora_counts <- c()
Snod_breadth_counts <- c()

for (SRR in colnames(pandora_output)) {
  Gilli_pandora_counts <- c(Gilli_pandora_counts, length(which(pandora_output[intersecting_genes_Gilli, SRR] == 1)))
  Gilli_breadth_counts <- c(Gilli_breadth_counts, length(which(breadth_out_0.5[intersecting_genes_Gilli, SRR] == 1)))
  
  Snod_pandora_counts <- c(Snod_pandora_counts, length(which(pandora_output[intersecting_genes_Snod, SRR] == 1)))
  Snod_breadth_counts <- c(Snod_breadth_counts, length(which(breadth_out_0.5[intersecting_genes_Snod, SRR] == 1)))
}

par(mfrow=c(2,1))
plot(Gilli_pandora_counts, Gilli_breadth_counts, main = "Gilliamella", xlab = "# pandora positive calls", ylab = "# breadth (50%) positive calls",
     xlim = c(0, 2250), ylim = c(0, 2250))
abline(b = 1, a = 0, lty = 2, lwd = 2)

plot(Snod_pandora_counts, Snod_breadth_counts, main = "Snodgrassella", xlab = "# pandora positive calls", ylab = "# breadth (50%) positive calls",
xlim = c(0, 1500), ylim = c(0, 1500))
abline(b = 1, a = 0, lty = 2, lwd = 2)
par(mfrow=c(1,1))



### For one sample (SRR7287243) compare depth and coverage of genes that intersect vs those that don't

Gilli_mean_depth <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_mean_depth_per_site.rds")

SRR <- "SRR7287243"

pandora_pos_genes <- intersecting_genes_Gilli[which(pandora_output[intersecting_genes_Gilli, SRR] == 1)]
breadth0.5_pos_genes <- intersecting_genes_Gilli[which(breadth_out_0.5[intersecting_genes_Gilli, SRR] == 1)]

pandora_only_genes <- pandora_pos_genes[which(! pandora_pos_genes %in% breadth0.5_pos_genes)]
breadth0.5_only_genes <- breadth0.5_pos_genes[which(! breadth0.5_pos_genes %in% pandora_pos_genes)]
both_genes <- pandora_pos_genes[which(pandora_pos_genes %in% breadth0.5_pos_genes)]

par(mfrow = c(2, 2))
boxplot(Gilli_mean_depth[pandora_only_genes, SRR], Gilli_mean_depth[breadth0.5_only_genes, SRR], Gilli_mean_depth[both_genes, SRR],
        main = "Gilli (sample SRR7287243)", names = c("pandora only", "breadth0.5 only", "both"), ylab = "Mean depth per site")

boxplot(breadth_out[pandora_only_genes, SRR], breadth_out[breadth0.5_only_genes, SRR], breadth_out[both_genes, SRR],
        main = "Gilli (sample SRR7287243)", names = c("pandora only", "breadth0.5 only", "both"), ylab = "Breadth of coverage")


Snod_mean_depth <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Snodgrassella_mean_depth_per_site.rds")

pandora_pos_genes <- intersecting_genes_Snod[which(pandora_output[intersecting_genes_Snod, SRR] == 1)]
breadth0.5_pos_genes <- intersecting_genes_Snod[which(breadth_out_0.5[intersecting_genes_Snod, SRR] == 1)]

pandora_only_genes <- pandora_pos_genes[which(! pandora_pos_genes %in% breadth0.5_pos_genes)]
breadth0.5_only_genes <- breadth0.5_pos_genes[which(! breadth0.5_pos_genes %in% pandora_pos_genes)]
both_genes <- pandora_pos_genes[which(pandora_pos_genes %in% breadth0.5_pos_genes)]


boxplot(Snod_mean_depth[pandora_only_genes, SRR], Snod_mean_depth[breadth0.5_only_genes, SRR], Snod_mean_depth[both_genes, SRR],
        main = "Snod (sample SRR7287243)", names = c("pandora only", "breadth0.5 only", "both"), ylab = "Mean depth per site")

boxplot(breadth_out[pandora_only_genes, SRR], breadth_out[breadth0.5_only_genes, SRR], breadth_out[both_genes, SRR],
        main = "Snod (sample SRR7287243)", names = c("pandora only", "breadth0.5 only", "both"), ylab = "Breadth of coverage")
par(mfrow = c(1, 1))



### Check which method better identifies core genes as present

### Gilliamella

Gilliamella_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Gilliamella_gene_presence_absence_roary-formatted.csv.gz"

Gilliamella_panaroo_out <- read.table(Gilliamella_path_to_panaroo,
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilli", rownames(Gilliamella_panaroo_out), sep = "_")

Gilliamella_panaroo_core_genes <- rownames(Gilliamella_panaroo_out)[which(Gilliamella_panaroo_out$No..isolates >= 99)]
prop_core_gene <- colSums(pandora_output[Gilliamella_panaroo_core_genes, ]) / colSums(pandora_output[grep("Gilli", rownames(Gilliamella_panaroo_out)), ])
colSums(pandora_output[grep("Gilli", rownames(Gilliamella_panaroo_out)), ])

breadth_core_genes_present <- Gilliamella_panaroo_core_genes[which(Gilliamella_panaroo_core_genes %in% rownames(breadth_out_0.5))]
colSums(breadth_out_0.5[breadth_core_genes_present, ])

par(mfrow = c(2, 1))
plot(colSums(pandora_output[grep("Gilli", rownames(Gilliamella_panaroo_out)), ]),
     colSums(pandora_output[Gilliamella_panaroo_core_genes, ]),
     ylim = c(0, 700),
     xlim = c(0, 9000),
     xlab = "Total Gilliamella genes",
     ylab = "Total core genes identified",
     main = "pandora")
abline(a = length(Gilliamella_panaroo_core_genes), b = 0, lwd = 1, lty = 2)


plot(colSums(breadth_out_0.5[grep("Gilli", rownames(breadth_out_0.5)), ]),
     colSums(breadth_out_0.5[breadth_core_genes_present, ]),
     ylim = c(0, 700),
     xlim = c(0, 9000),
     xlab = "Total Gilliamella genes",
     ylab = "Total core genes identified",
     main = "breadth (0.5)")
abline(a = length(Gilliamella_panaroo_core_genes), b = 0, lwd = 1, lty = 2)

par(mfrow = c(1, 1))


### Snodgrassella

Snodgrassella_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Snodgrassella_gene_presence_absence_roary-formatted.csv.gz"

Snodgrassella_panaroo_out <- read.table(Snodgrassella_path_to_panaroo,
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snod", rownames(Snodgrassella_panaroo_out), sep = "_")

Snodgrassella_panaroo_core_genes <- rownames(Snodgrassella_panaroo_out)[which(Snodgrassella_panaroo_out$No..isolates >= 53)]
prop_core_gene <- colSums(pandora_output[Snodgrassella_panaroo_core_genes, ]) / colSums(pandora_output[grep("Snod", rownames(Snodgrassella_panaroo_out)), ])
colSums(pandora_output[grep("Snod", rownames(Snodgrassella_panaroo_out)), ])

breadth_core_genes_present <- Snodgrassella_panaroo_core_genes[which(Snodgrassella_panaroo_core_genes %in% rownames(breadth_out_0.5))]
colSums(breadth_out_0.5[breadth_core_genes_present, ])

par(mfrow = c(2, 1))
plot(colSums(pandora_output[grep("Snod", rownames(Snodgrassella_panaroo_out)), ]),
     colSums(pandora_output[Snodgrassella_panaroo_core_genes, ]),
     ylim = c(0, 800),
     xlim = c(0, 4000),
     xlab = "Total Snodgrassella genes",
     ylab = "Total core genes identified",
     main = "pandora")
abline(a = length(Snodgrassella_panaroo_core_genes), b = 0, lwd = 1, lty = 2)


plot(colSums(breadth_out_0.5[grep("Snod", rownames(breadth_out_0.5)), ]),
     colSums(breadth_out_0.5[breadth_core_genes_present, ]),
     ylim = c(0, 800),
     xlim = c(0, 4000),
     xlab = "Total Snodgrassella genes",
     ylab = "Total core genes identified",
     main = "breadth (0.5)")
abline(a = length(Snodgrassella_panaroo_core_genes), b = 0, lwd = 1, lty = 2)

par(mfrow = c(1, 1))




### Bifidobacterium

Bifidobacterium_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/panaroo_pangenomes/complete_raw_output/Bifidobacterium/gene_presence_absence_roary.csv.gz"

Bifidobacterium_panaroo_out <- read.table(Bifidobacterium_path_to_panaroo,
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Bifidobacterium_panaroo_out) <- paste("Bifidobacterium", rownames(Bifidobacterium_panaroo_out), sep = "_")

Bifidobacterium_panaroo_core_genes <- rownames(Bifidobacterium_panaroo_out)[which(Bifidobacterium_panaroo_out$No..isolates >= 14)]
prop_core_gene <- colSums(pandora_output[Bifidobacterium_panaroo_core_genes, ]) / colSums(pandora_output[grep("Bifidobacterium", rownames(Bifidobacterium_panaroo_out)), ])
colSums(pandora_output[grep("Bifidobacterium", rownames(Bifidobacterium_panaroo_out)), ])

plot(colSums(pandora_output[grep("Bifidobacterium", rownames(Bifidobacterium_panaroo_out)), ]),
     colSums(pandora_output[Bifidobacterium_panaroo_core_genes, ]),
     ylim = c(0, 1200),
     xlim = c(0, 1500),
     xlab = "Total Bifidobacterium genes",
     ylab = "Total core genes identified",
     main = "pandora - Bifidobacterium")
abline(a = length(Bifidobacterium_panaroo_core_genes), b = 0, lwd = 1, lty = 2)



### Firm4

Firm4_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/panaroo_pangenomes/complete_raw_output/Firm4//gene_presence_absence_roary.csv.gz"

Firm4_panaroo_out <- read.table(Firm4_path_to_panaroo,
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Firm4_panaroo_out) <- paste("Firm4", rownames(Firm4_panaroo_out), sep = "_")

Firm4_panaroo_core_genes <- rownames(Firm4_panaroo_out)[which(Firm4_panaroo_out$No..isolates >= 6)]
prop_core_gene <- colSums(pandora_output[Firm4_panaroo_core_genes, ]) / colSums(pandora_output[grep("Firm4", rownames(Firm4_panaroo_out)), ])
colSums(pandora_output[grep("Firm4", rownames(Firm4_panaroo_out)), ])

plot(colSums(pandora_output[grep("Firm4", rownames(Firm4_panaroo_out)), ]),
     colSums(pandora_output[Firm4_panaroo_core_genes, ]),
     ylim = c(0, 50),
     xlim = c(0, 4000),
     xlab = "Total Firm4 genes",
     ylab = "Total core genes identified",
     main = "pandora - Firm4")
abline(a = length(Firm4_panaroo_core_genes), b = 0, lwd = 1, lty = 2)


### Firm5

Firm5_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/panaroo_pangenomes/complete_raw_output/Firm5//gene_presence_absence_roary.csv.gz"

Firm5_panaroo_out <- read.table(Firm5_path_to_panaroo,
                                header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Firm5_panaroo_out) <- paste("Firm5", rownames(Firm5_panaroo_out), sep = "_")

Firm5_panaroo_core_genes <- rownames(Firm5_panaroo_out)[which(Firm5_panaroo_out$No..isolates >= 25)]
prop_core_gene <- colSums(pandora_output[Firm5_panaroo_core_genes, ]) / colSums(pandora_output[grep("Firm5", rownames(Firm5_panaroo_out)), ])
colSums(pandora_output[grep("Firm5", rownames(Firm5_panaroo_out)), ])

plot(colSums(pandora_output[grep("Firm5", rownames(Firm5_panaroo_out)), ]),
     colSums(pandora_output[Firm5_panaroo_core_genes, ]),
     ylim = c(0, 300),
     xlim = c(0, 4500),
     xlab = "Total Firm5 genes",
     ylab = "Total core genes identified",
     main = "pandora - Firm5")
abline(a = length(Firm5_panaroo_core_genes), b = 0, lwd = 1, lty = 2)
