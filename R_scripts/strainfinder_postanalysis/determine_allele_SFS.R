rm(list = ls(all.names = TRUE))

# Determine SFS for all haplotyped alleles
# Save in convenient format.

all_genes_abun <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/all_genes_OTU.rds")

raw_all_gene_sfs <- list()

for (sp in names(all_genes_abun)) {

  raw_sp_gene_sfs <- list()
  
  for (g in names(all_genes_abun[[sp]])) {
      
     allele_names <- paste(g, as.character(1:ncol(all_genes_abun[[sp]][[g]])), sep = "|||")
     
     allele_occurences <- as.integer(colSums(all_genes_abun[[sp]][[g]] > 0))
    
     raw_sp_gene_sfs[[g]] <- data.frame(Species = sp,
                                        Gene = g,
                                        Allele = allele_names,
                                        Sample_occurrences = allele_occurences)
     
  }
  
  raw_all_gene_sfs[[sp]] <- do.call(rbind, raw_sp_gene_sfs)

}

all_gene_sfs <- do.call(rbind, raw_all_gene_sfs)

write.table(x = all_gene_sfs,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/allele_sample_occurrences.tsv",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)

