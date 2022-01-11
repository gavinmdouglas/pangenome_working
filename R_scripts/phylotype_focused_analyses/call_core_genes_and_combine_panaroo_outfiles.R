### Call core genes for each phylotype and make a single table with all genes mapped to.

rm(list = ls(all.names = TRUE))

core_genes <- list()

Bifidobacterium_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Bifidobacterium/gene_presence_absence_roary.csv",
                                      header = TRUE, sep = ",", comment.char = "", quote = "")
rownames(Bifidobacterium_panaroo) <- paste("Bifidobacterium", Bifidobacterium_panaroo$Gene, sep = "_")
core_genes[["Bifidobacterium"]] <- rownames(Bifidobacterium_panaroo)[which(Bifidobacterium_panaroo$No..isolates >= 14)]


Firm4_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Firm4/gene_presence_absence_roary.csv.gz",
                            header = TRUE, sep = ",", comment.char = "", quote = "")
rownames(Firm4_panaroo) <- paste("Firm4", Firm4_panaroo$Gene, sep = "_")
core_genes[["Firm4"]] <- rownames(Firm4_panaroo)[which(Firm4_panaroo$No..isolates >= 5)]


Firm5_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Firm5/gene_presence_absence_roary.csv",
                            header = TRUE, sep = ",", comment.char = "", quote = "")
rownames(Firm5_panaroo) <- paste("Firm5", Firm5_panaroo$Gene, sep = "_")
core_genes[["Firm5"]] <- rownames(Firm5_panaroo)[which(Firm5_panaroo$No..isolates >= 25)]


Gilliamella_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Gilliamella/gene_presence_absence_roary.csv",
                                  header = TRUE, sep = ",", comment.char = "", quote = "")
rownames(Gilliamella_panaroo) <- paste("Gilliamella", Gilliamella_panaroo$Gene, sep = "_")
core_genes[["Gilliamella"]] <- rownames(Gilliamella_panaroo)[which(Gilliamella_panaroo$No..isolates >= 98)]


Snodgrassella_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Snodgrassella/gene_presence_absence_roary.csv",
                                    header = TRUE, sep = ",", comment.char = "", quote = "")
rownames(Snodgrassella_panaroo) <- paste("Snodgrassella", Snodgrassella_panaroo$Gene, sep = "_")
core_genes[["Snodgrassella"]] <- rownames(Snodgrassella_panaroo)[which(Snodgrassella_panaroo$No..isolates >= 51)]

all_core_genes <- c(core_genes$Bifidobacterium, core_genes$Firm4, core_genes$Firm5, core_genes$Gilliamella, core_genes$Snodgrassella)


# Combine all panaroo tables into single master dataframe
all_panaroo <- rbind(Bifidobacterium_panaroo[, c(1:14)], Firm4_panaroo[, c(1:14)])
all_panaroo <- rbind(all_panaroo, Firm5_panaroo[, c(1:14)])
all_panaroo <- rbind(all_panaroo, Gilliamella_panaroo[, c(1:14)])
all_panaroo <- rbind(all_panaroo, Snodgrassella_panaroo[, c(1:14)])

panaroo_summary <- list(core_genes = all_core_genes, panaroo_all_phylotypes = all_panaroo)

saveRDS(object = panaroo_summary, file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/combined_panaroo_and_core_genes.rds")

