Apilactobacillus_apinorum_ids <- read.table("/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes_prokka/Apilactobacillus_apinorum/GCA_001281175.1/GCA_001281175.1.ids.txt",
                                            header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1
sp_panaroo <- data.frame(matrix(NA, nrow = length(Apilactobacillus_apinorum_ids), ncol = 15))
colnames(sp_panaroo) <- c("Gene", "Non.unique.Gene.name", "Annotation", "No..isolates",
                          "No..sequences", "Avg.sequences.per.isolate", "Genome.Fragment",
                          "Order.within.Fragment", "Accessory.Fragment", "Accessory.Order.with.Fragment",
                          "QC","Min.group.size.nuc", "Max.group.size.nuc", "Avg.group.size.nuc", "GCA_001281175.1")

sp_panaroo$GCA_001281175.1 <- Apilactobacillus_apinorum_ids
sp_panaroo$Gene <- Apilactobacillus_apinorum_ids

write.table(x = sp_panaroo,
            file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Apilactobacillus_apinorum/gene_presence_absence_roary.csv",
            sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

sp_panaroo_simple <- sp_panaroo[, c("Gene", "Non.unique.Gene.name", "Annotation", "GCA_001281175.1")]
write.table(x = sp_panaroo_simple,
            file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Apilactobacillus_apinorum/gene_presence_absence.csv",
            sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

