rm(list = ls(all.names = TRUE))

# Identify strain that produces the highest Jaccard for each allele (based on # intersecting samples).

library(parallel)

strain_abun <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/strain_abun.rds")

all_genes_abun <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/all_genes_OTU.rds")

presence_matrix <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.matrix",
                              header = TRUE, sep = "\t", row.names = 1)

species_presence <- readRDS(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/species_presence_50core_or_10percent.rds")


all_genes <- unlist(sapply(all_genes_abun, names))


identify_strain_allele_combos_w_highest_jaccard <- function(gene_name) {

   gene_split <- strsplit(gene_name, "_")[[1]]
   
   sp <- paste(gene_split[1], gene_split[2], sep = "_")
   
   if (sp == "Bifidobacterium_coryneforme") {
      sp <- "Bifidobacterium_coryneforme_indicum"
   }
   
   allele_abun <- all_genes_abun[[sp]][[gene_name]]
   sp_strains <- strain_abun[[sp]]
   
   colnames(allele_abun) <- paste("haplotype", as.character(1:ncol(allele_abun)), sep = "_")
   
   out_df <- data.frame(matrix(NA, nrow = ncol(allele_abun), ncol = 6))
   colnames(out_df) <- c("Species", "Gene", "Allele", "Jaccard", "Num_best_matches", "Best_matches")
   out_df$Species <- sp
   out_df$Gene <- gene_name
   out_df$Allele <- colnames(allele_abun)
   out_df$Jaccard <- 0
   out_df$Num_best_matches <- 0
   out_df$Best_matches <- NA

   intersecting_samples <- rownames(allele_abun)[which(rownames(allele_abun) %in% rownames(sp_strains))]
   
   if (length(intersecting_samples) == 0) { return(out_df) }
   
   allele_abun <- allele_abun[intersecting_samples, , drop = FALSE]
   
   for (allele_i in 1:ncol(allele_abun)) {
    
      allele_samples <- rownames(allele_abun)[which(allele_abun[, allele_i] > 0)]
      
      all_strain_jaccard <- as.numeric()
      
      for (strain_i in 1:ncol(sp_strains)) {
         strain_samples <- rownames(sp_strains)[which(sp_strains[, strain_i] > 0)]
         
         num_intersect_samples <- length(which(allele_samples %in% strain_samples))
         
         num_union_samples <- length(unique(c(allele_samples, strain_samples)))
         
         all_strain_jaccard <- c(all_strain_jaccard, num_intersect_samples / num_union_samples)
      }

      max_jaccard <- max(all_strain_jaccard)
      
      out_df[allele_i, "Jaccard"] <- max_jaccard
      out_df[allele_i, "Num_best_matches"] <- length(which(all_strain_jaccard == max_jaccard))
      out_df[allele_i, "Best_matches"] <- paste(colnames(sp_strains)[which(all_strain_jaccard == max_jaccard)], collapse = ";")
      
   }

   return(out_df)
   
}


raw_jaccard_comparisons <- mclapply(all_genes, identify_strain_allele_combos_w_highest_jaccard, mc.cores = 50)

jaccard_comparisons <- do.call(rbind, raw_jaccard_comparisons)

saveRDS(object = jaccard_comparisons,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/strain_vs_gene_abun_best_jaccard.rds")




### TESTING

library(ggplot2)
library(ggbeeswarm)

jaccard_comparisons_mean <- aggregate(Jaccard ~ Species + Gene,
                                      data = jaccard_comparisons,
                                      FUN = mean)

jaccard_comparisons_mean_by_species <- aggregate(Jaccard ~ Species,
                                                 data = jaccard_comparisons_mean,
                                                 FUN = mean)

jaccard_comparisons_mean_by_species <- jaccard_comparisons_mean_by_species[order(jaccard_comparisons_mean_by_species$Jaccard), ]

jaccard_comparisons_mean$Species <- factor(jaccard_comparisons_mean$Species, levels = jaccard_comparisons_mean_by_species$Species)

ggplot(data = jaccard_comparisons_mean, aes(x = Jaccard, y = Species)) +
   geom_quasirandom(groupOnX = FALSE) +
   geom_boxplot(alpha = 0.6, fill = "grey", outlier.shape = NA) +
   theme_bw() +
   xlab("Highest possible Jaccard similarity between allele and strain presence/absence\n(averaged across all alleles per gene)")


# Get quick stats of # strains, mean # alleles, and the num core, and num accessory genes per species.

species_info <- data.frame(num_samples_w_strains = sapply(strain_abun, nrow),
                           num_strains = sapply(strain_abun, ncol))

species_info$mean_num_alleles <- sapply(names(all_genes_abun), function(x) { mean(sapply(all_genes_abun[[x]], ncol)) })

species_info$total_genes <- NA
species_info$noncore_genes <- NA

presence_matrix_noncore <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.noncore.matrix",
                                      header = TRUE, sep = "\t", row.names = 1)

for (sp in rownames(species_info)) {
   species_info[sp, "total_genes"] <- length(grep(sp, rownames(presence_matrix)))
   species_info[sp, "noncore_genes"] <- length(grep(sp, rownames(presence_matrix_noncore)))
}

species_info$core_genes <- species_info$total_genes - species_info$noncore_genes
