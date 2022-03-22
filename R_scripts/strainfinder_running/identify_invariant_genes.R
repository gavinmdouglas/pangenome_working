# Identify genes that are called as present, but which lack mutations called by pandora.

rm(list = ls(all.names = TRUE))

presence_matrix <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.matrix",
                               header = TRUE, sep = "\t", row.names = 1)

sf_prepped_genes <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/all_genes_prepped.txt",
                               stringsAsFactors = FALSE, header = FALSE)$V1

presence_matrix_num_samples <- rowSums(presence_matrix > 0)

sf_prepped_genes <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/all_genes_prepped.txt",
                               stringsAsFactors = FALSE, header = FALSE)$V1

boxplot(presence_matrix_num_samples, presence_matrix_num_samples[-which(names(presence_matrix_num_samples) %in% sf_prepped_genes)],
        names = c("All genes", "Genes without variants"), ylab = "Number of samples")

invariant_genes <- names(presence_matrix_num_samples[-which(names(presence_matrix_num_samples) %in% sf_prepped_genes)])

write.table(x = invariant_genes, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/invariant_genes.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)
