rm(list = ls(all.names = TRUE))

# Get measure of how well gene and strain haplotypes match up in order to see if they correlate with gene tree DTL metrics.
# Based on RVadj metric comparing matrices, which is an extension of Pearson correlation.

library(MatrixCorrelation)
library(parallel)

set.seed(20220321)

strain_abun <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/strain_abun.rds")

all_genes_abun <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/all_genes_OTU.rds")

presence_matrix <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.matrix",
                              header = TRUE, sep = "\t", row.names = 1)

species_presence <- readRDS(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/species_presence_50core_or_10percent.rds")


all_genes <- unlist(sapply(all_genes_abun, names))



create_RVadj_row <- function(gene) {
 
   gene_split <- strsplit(gene, "_")[[1]]
  
   sp <- paste(gene_split[1], gene_split[2], sep = "_")
   
   if (sp == "Bifidobacterium_coryneforme") {
     sp <- "Bifidobacterium_coryneforme_indicum"
   }
  
   species_samples <- colnames(species_presence)[which(species_presence[sp, ] > 0)]
   
   species_samples_w_strain <- rownames(strain_abun[[sp]])
   
   species_samples_wo_strain <- species_samples[which(! species_samples %in% species_samples_w_strain)]
   
   gene_samples <- colnames(presence_matrix)[which(presence_matrix[gene, ] > 0)]
   
   gene_samples_w_haplotype <- rownames(all_genes_abun[[sp]][[gene]])
   
   gene_samples_wo_haplotype <- gene_samples[which(! gene_samples %in% gene_samples_w_haplotype)]
   
   intersecting_samples <- species_samples_w_strain[which(species_samples_w_strain %in% gene_samples_w_haplotype)]
   
   samples2exclude <- c(species_samples_wo_strain, gene_samples_wo_haplotype)

   strain_abun_subset <- strain_abun[[sp]][intersecting_samples, , drop = FALSE]
   gene_abun_subset <- all_genes_abun[[sp]][[gene]][intersecting_samples, , drop = FALSE]
  
   sample_concordance <- data.frame(matrix(NA, nrow = 1, ncol = 9))
   colnames(sample_concordance) <- c("species", "gene", "num_strain_samples", "num_gene_samples", "num_intersecting",
                                     "num_excluded_samples", "RVadj", "mean_permuted_RVadj", "permuted_RVadj_p")
   
   sample_concordance[1, c("species", "gene")] <- c(sp, gene)
   sample_concordance[1, c("num_strain_samples",
                           "num_gene_samples",
                           "num_intersecting",
                           "num_excluded_samples")] <- c(length(species_samples_w_strain),
                                                         length(gene_samples_w_haplotype),
                                                         length(intersecting_samples),
                                                         length(samples2exclude))
   
   if (length(intersecting_samples) <= 2) { 
     return(sample_concordance)
   }
   
   sample_concordance[1, "RVadj"] <- MatrixCorrelation::RVadj(as.matrix(strain_abun_subset),
                                                              as.matrix(gene_abun_subset))

   if (length(intersecting_samples) >= 20) {
     rand_RVadj <- as.numeric()
     
     for (i in 1:999) {
       
       gene_abun_subset_rand <- gene_abun_subset
       
       for (j in 1:nrow(gene_abun_subset_rand)) {
         gene_abun_subset_rand[j, ] <- sample(gene_abun_subset[j, , drop = FALSE])
       }
       
       rand_RVadj <- c(rand_RVadj, MatrixCorrelation::RVadj(as.matrix(strain_abun_subset),
                                                            as.matrix(gene_abun_subset_rand)))
      
     }
     
     sample_concordance[1, "mean_permuted_RVadj"] <- mean(rand_RVadj)
     
     sample_concordance[1, "permuted_RVadj_p"] <- (length(which(rand_RVadj >= sample_concordance[1, "RVadj"])) + 1) / 1000

   }
   
   return(sample_concordance)
   
}


RVadj_rows <- mclapply(all_genes, create_RVadj_row, mc.cores = 50)

RVadj_out <- do.call(rbind, RVadj_rows)

saveRDS(object = RVadj_out,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/strain_vs_gene_abun_RVadj.rds")

