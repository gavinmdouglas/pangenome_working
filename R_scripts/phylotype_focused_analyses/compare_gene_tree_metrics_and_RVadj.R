rm(list = ls(all.names = TRUE))

gene_strain_abun_RVadj <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/phylotype_prepped_files/all_individual_genes/output/gene_strain_abun_RVadj.rds")

gene_tree_metrics <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/phylotype/gene_tree_analysis/ranger_dtl_output/rangerdtl_mean_event_counts.tsv",
                                header = FALSE, sep = "\t", row.names = 1)

colnames(gene_tree_metrics) <- c("dup", "transfer", "loss")

gene_tree_metrics$summed <- rowSums(gene_tree_metrics)

intersecting_genes <- rownames(gene_strain_abun_RVadj)[which(rownames(gene_strain_abun_RVadj) %in% rownames(gene_tree_metrics))]

cor.test(gene_tree_metrics[intersecting_genes, "transfer"],
         gene_strain_abun_RVadj[intersecting_genes, "RVadj"],
         method="spearman")
plot(gene_tree_metrics[intersecting_genes, "transfer"],
     gene_strain_abun_RVadj[intersecting_genes, "RVadj"],
     xlab = "transfer DTL",
     ylab = "RVadj",
     ylim = c(-0.3, 1),
     main = "All phylotypes")


par(mfrow = c(2, 2))
intersecting_genes_Bifidobacterium <- grep("Bifidobacterium", intersecting_genes, value=TRUE)
cor.test(gene_tree_metrics[intersecting_genes_Bifidobacterium, "transfer"],
         gene_strain_abun_RVadj[intersecting_genes_Bifidobacterium, "RVadj"],
         method="spearman")
plot(gene_tree_metrics[intersecting_genes_Bifidobacterium, "transfer"],
     gene_strain_abun_RVadj[intersecting_genes_Bifidobacterium, "RVadj"],
     xlab = "transfer DTL",
     ylab = "RVadj",
     ylim = c(-0.3, 1),
     main = "Bifidobacterium")


intersecting_genes_Firm5 <- grep("Firm5", intersecting_genes, value=TRUE)
cor.test(gene_tree_metrics[intersecting_genes_Firm5, "transfer"],
         gene_strain_abun_RVadj[intersecting_genes_Firm5, "RVadj"],
         method="spearman")
plot(gene_tree_metrics[intersecting_genes_Firm5, "transfer"],
     gene_strain_abun_RVadj[intersecting_genes_Firm5, "RVadj"],
     xlab = "transfer DTL",
     ylab = "RVadj",
     ylim = c(-0.3, 1),
     main = "Firm5")



intersecting_genes_Gilliamella <- grep("Gilliamella", intersecting_genes, value=TRUE)
cor.test(gene_tree_metrics[intersecting_genes_Gilliamella, "transfer"],
         gene_strain_abun_RVadj[intersecting_genes_Gilliamella, "RVadj"],
         method="spearman")
plot(gene_tree_metrics[intersecting_genes_Gilliamella, "transfer"],
     gene_strain_abun_RVadj[intersecting_genes_Gilliamella, "RVadj"],
     xlab = "transfer DTL",
     ylab = "RVadj",
     ylim = c(-0.3, 1),
     main = "Gilliamella")



intersecting_genes_Snodgrassella <- grep("Snodgrassella", intersecting_genes, value=TRUE)
cor.test(gene_tree_metrics[intersecting_genes_Snodgrassella, "transfer"],
         gene_strain_abun_RVadj[intersecting_genes_Snodgrassella, "RVadj"],
         method="spearman")
plot(gene_tree_metrics[intersecting_genes_Snodgrassella, "transfer"],
     gene_strain_abun_RVadj[intersecting_genes_Snodgrassella, "RVadj"],
     xlab = "transfer DTL",
     ylab = "RVadj",
     ylim = c(-0.3, 1),
     main = "Snodgrassella")

par(mfrow = c(1, 1))