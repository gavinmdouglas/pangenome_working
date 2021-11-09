### Summary counts and functional annotations of identified clusters

rm(list = ls(all.names = TRUE))

library(parallel)

summarize_panaroo_annotations <- function(gene_set, panaroo_table) {
  
  all_annot <- panaroo_table[gene_set, "Annotation"]
  
  if (length(which(duplicated(all_annot))) > 0) {
    all_annot_unique <- all_annot[-which(duplicated(all_annot))]
    
    for (i in 1:length(all_annot_unique)) {
     
      num_matches <- length(which(all_annot == all_annot_unique[i]))
      
      if (num_matches > 1) {
        all_annot_unique[i] <- paste(all_annot_unique[i], " (x ", as.character(num_matches), ")", sep = "")
      }
       
    }
    
  } else {
    all_annot_unique <- all_annot 
  }
  
  return(paste(sort(all_annot_unique), collapse = " ||| "))
  
}

Gilliamella_panaroo <- read.table(file = "projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Gilliamella_gene_presence_absence_roary-formatted.csv.gz",
                                  header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE, quote = "", comment.char = "")
rownames(Gilliamella_panaroo) <- paste("Gilli", rownames(Gilliamella_panaroo), sep = "_")

Gilliamella_noncore_same_strain_cooccur_clusters <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Gilliamella_noncore_same_strain_cooccur_clusters.rds")

hist(sapply(Gilliamella_noncore_same_strain_cooccur_clusters, length), breaks = 100, main = "Gilliamella", xlab = "Number of genes per cluster")

Gilliamella_cluster_annot <- lapply(Gilliamella_noncore_same_strain_cooccur_clusters, summarize_panaroo_annotations, panaroo_table = Gilliamella_panaroo)


# Get enriched annotations for all clusters with < 50 members and then for each individual cluster with > 50 members.

fishers_exact_sig_and_background <- function(focal_func, sig_funcs, background_funcs) {
 
  focal_func <- "HTH-type transcriptional regulator LutR"
  sig_funcs <- Gilliamella_noncore_same_strain_cooccur_clusters[[110]]
  background_funcs <- Gilliamella_panaroo$Annotation
  
  focal_matrix <- matrix(c(length(which(sig_funcs == focal_func)),
                           length(which(sig_funcs != focal_func)),
                           length(which(background_funcs == focal_func)),
                           length(which(background_funcs != focal_func))),
                         nrow = 2, ncol = 2) 
  
}
 

enriched_func_parallel <- function(sig_funcs, background_funcs, num_cores = 25) {

  unique_background_funcs <- background_funcs[-which(duplicated(background_funcs))]
  
  all_enrichment_tests <- mclapply(X = unique_background_funcs, FUN = fishers_exact_sig_and_background, sig_funcs = sig_funcs, background_funcs = background_funcs, mc.cores = num_cores)

}


Snodgrassella_panaroo <- read.table(file = "projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Snodgrassella_gene_presence_absence_roary-formatted.csv.gz",
                                  header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE, quote = "", comment.char = "")
rownames(Snodgrassella_panaroo) <- paste("Snod", rownames(Snodgrassella_panaroo), sep = "_")

Snodgrassella_noncore_same_strain_cooccur_clusters <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_noncore_genome_and_breadth_cooccur/Snodgrassella_noncore_same_strain_cooccur_clusters.rds")

hist(sapply(Snodgrassella_noncore_same_strain_cooccur_clusters, length), breaks = 100, main = "Snodgrassella", xlab = "Number of genes per cluster")
