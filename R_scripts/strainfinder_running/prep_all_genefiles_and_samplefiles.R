### Create genefiles and samplefiles for each individual gene that is present in at least one sample.

rm(list = ls(all.names = TRUE))

gene_matrix <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.matrix",
                          header = TRUE, sep = "\t", row.names = 1)

setwd("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/starting_files")

for (gene in rownames(gene_matrix)) {

  genefile <- paste(gene, "_genes.txt", sep = "") 
  write.table(gene, file = genefile, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  samplefile <- paste(gene, "_samples.txt", sep = "")
  samples_with_gene <- colnames(gene_matrix)[which(gene_matrix[gene, ] > 0)]
  write.table(samples_with_gene, file = samplefile, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
}

write.table(x = rownames(gene_matrix), file = "../all_genes.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
