### Remove samples with insufficient coverage per genus from the input file for pandora.

rm(list = ls(all.names = TRUE))

pandora_output <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_multisample.matrix",
                             header = TRUE, sep = "\t", row.names = 1)
rownames(pandora_output) <- gsub("\\.fa$", "", rownames(pandora_output))

### Get core genes:
Bifidobacterium_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/panaroo_pangenomes/complete_raw_output/Bifidobacterium/gene_presence_absence_roary.csv.gz"

Bifidobacterium_panaroo_out <- read.table(Bifidobacterium_path_to_panaroo,
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Bifidobacterium_panaroo_out) <- paste("Bifidobacterium", rownames(Bifidobacterium_panaroo_out), sep = "_")

Bifidobacterium_panaroo_core_genes <- rownames(Bifidobacterium_panaroo_out)[which(Bifidobacterium_panaroo_out$No..isolates >= 14)]



Firm4_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/panaroo_pangenomes/complete_raw_output/Firm4//gene_presence_absence_roary.csv.gz"

Firm4_panaroo_out <- read.table(Firm4_path_to_panaroo,
                                header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Firm4_panaroo_out) <- paste("Firm4", rownames(Firm4_panaroo_out), sep = "_")

Firm4_panaroo_core_genes <- rownames(Firm4_panaroo_out)[which(Firm4_panaroo_out$No..isolates >= 6)]



Firm5_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/panaroo_pangenomes/complete_raw_output/Firm5//gene_presence_absence_roary.csv.gz"

Firm5_panaroo_out <- read.table(Firm5_path_to_panaroo,
                                header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Firm5_panaroo_out) <- paste("Firm5", rownames(Firm5_panaroo_out), sep = "_")

Firm5_panaroo_core_genes <- rownames(Firm5_panaroo_out)[which(Firm5_panaroo_out$No..isolates >= 25)]



Gilliamella_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Gilliamella_gene_presence_absence_roary-formatted.csv.gz"

Gilliamella_panaroo_out <- read.table(Gilliamella_path_to_panaroo,
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilli", rownames(Gilliamella_panaroo_out), sep = "_")

Gilliamella_panaroo_core_genes <- rownames(Gilliamella_panaroo_out)[which(Gilliamella_panaroo_out$No..isolates >= 99)]


Snodgrassella_path_to_panaroo <- "/data1/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Snodgrassella_gene_presence_absence_roary-formatted.csv.gz"

Snodgrassella_panaroo_out <- read.table(Snodgrassella_path_to_panaroo,
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snod", rownames(Snodgrassella_panaroo_out), sep = "_")

Snodgrassella_panaroo_core_genes <- rownames(Snodgrassella_panaroo_out)[which(Snodgrassella_panaroo_out$No..isolates >= 53)]

core_genes <- list("Bifidobacterium" = Bifidobacterium_panaroo_core_genes,
                   "Firm4" = Firm4_panaroo_core_genes,
                   "Firm5" = Firm5_panaroo_core_genes,
                   "Gilli" = Gilliamella_panaroo_core_genes,
                   "Snod" = Snodgrassella_panaroo_core_genes)

# Table representing % of core genes called as present per phylotype per sample.
core_gene_coverage <- data.frame(matrix(NA, ncol = 5, nrow = 74))
colnames(core_gene_coverage) <- c("Bifidobacterium", "Firm4", "Firm5", "Gilli", "Snod")
rownames(core_gene_coverage) <- colnames(pandora_output)

for (phylotype in colnames(core_gene_coverage)) {
  for (samplename in rownames(core_gene_coverage)) {
    core_gene_coverage[samplename, phylotype] <- (length(which(pandora_output[core_genes[[phylotype]], samplename] == 1)) / length(core_genes[[phylotype]])) * 100
  }   
}

max_coverage_per_sample <- apply(core_gene_coverage, 1, max)

# Remove three obvious outliers:
plot(max_coverage_per_sample, colSums(pandora_output))

colnames(pandora_output)[which(colSums(pandora_output) < 2000)]

pandora_output_filt <- pandora_output[, -which(colSums(pandora_output) < 2000)]

write.table(x = pandora_output_filt,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_multisample.filt.matrix",
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

panaroo_core_genes <- c(Bifidobacterium_panaroo_core_genes, Firm4_panaroo_core_genes, Firm5_panaroo_core_genes,
                        Gilliamella_panaroo_core_genes, Snodgrassella_panaroo_core_genes)

pandora_output_filt_noncore <- pandora_output_filt[-which(rownames(pandora_output_filt) %in% panaroo_core_genes), ]

write.table(x = pandora_output_filt_noncore,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_multisample.filt.noncore.matrix",
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

