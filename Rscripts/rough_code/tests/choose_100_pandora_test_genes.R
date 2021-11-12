### Choose 100 test genes for pandora. Should be representative of different levels of prevalence and phylotypes.
### Write out gene file per gene as well as samples where the gene was called as present.

gene_matrix <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina/pandora_multisample.filt.matrix",
                          header = TRUE, sep = "\t", row.names = 1)


gene_prev <- rowSums(gene_matrix)

test_genes <- c()


rare_genes <- names(which((gene_prev > 1) & (gene_prev < 20)))

common_genes <- names(which(gene_prev > 20))


set.seed(141)
test_genes <- c(test_genes, sample(names(which(gene_prev == 1)), 20))
test_genes <- c(test_genes, sample(rare_genes, 40))
test_genes <- c(test_genes, sample(common_genes, 40))


setwd("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/starting_files/")

for (test_gene in test_genes) {
 
  genefile <- paste(test_gene, "_genefile.txt", sep = "") 
  write.table(test_gene, file = genefile, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  samplefile <- paste(test_gene, "_samples.txt", sep = "")
  samples_with_gene <- colnames(gene_matrix)[which(gene_matrix[test_gene, ] > 0)]
  write.table(samples_with_gene, file = samplefile, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
}

write.table(x = test_genes, file = "../all_test_genes.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)