rm(list = ls(all.names = TRUE))

# Get measure of how well gene and strain haplotypes match up in order to see if they correlate with gene tree DTL metrics.
# Based on RVadj metric comparing matrices, which is an extension of Pearson correlation.

library(MatrixCorrelation)

total_num_samples <- 71

core_haplotype_abun <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/phylotype_prepped_files/core_genes_output/strainfinder_pandora_struct_abun.rds")

all_gene_otu_tabs <- list.files(path = "//data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/phylotype_prepped_files/all_individual_genes/output/otu_tables",
                               pattern = "otu_table",
                               full.names = TRUE)

all_gene_otu_tabs_no_path <- list.files(path = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/phylotype_prepped_files/all_individual_genes/output/otu_tables",
                                pattern = "otu_table",
                                full.names = FALSE)

all_genes <- gsub("otu_table.", "", all_gene_otu_tabs_no_path)
all_genes <- gsub(".txt", "", all_genes)
all_genes <- gsub("\\..*$", "", all_genes)

sample_concordance <- data.frame(matrix(NA, nrow = length(all_genes), ncol = 6))
colnames(sample_concordance) <- c("phylotype", "gene", "num_gene_samples", "num_strain_samples", "num_intersecting", "RVadj")

rownames(sample_concordance) <- all_genes
sample_concordance$gene <- all_genes


for (i in 1:length(all_genes)) {

  # wc_cmd <- paste("wc -l", all_gene_samples[i])
  # if (system(wc_cmd)[1] == 0) { next }

  gene <- all_genes[i]
  
  otu_tab <- read.table(file = all_gene_otu_tabs[i], stringsAsFactors = FALSE, skip = 1)
  
  gene_samples_file <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/phylotype_prepped_files/all_individual_genes/prepped_input_w_new_struct_alleles/",
                             gene,
                             "_samples.txt", sep = "")
  
  gene_samples <- read.table(file = gene_samples_file, stringsAsFactors = FALSE)$V1
  

  phylotype <- strsplit(gene, split = "_")[[1]][1]
  phylotype_samples <- rownames(core_haplotype_abun[[phylotype]])
  
  intersecting_samples <- phylotype_samples[which(phylotype_samples %in% gene_samples)]
  
  sample_concordance[gene, "phylotype"] <- phylotype
  sample_concordance[gene, c("num_gene_samples", "num_strain_samples", "num_intersecting")] <- c(length(gene_samples),
                                                                                              length(phylotype_samples),
                                                                                              length(intersecting_samples))
  
  if (length(intersecting_samples) <= 2) { next }

  otu_tab_subset <- otu_tab[which(gene_samples %in% intersecting_samples), ]
  strain_tab_subset <- core_haplotype_abun[[phylotype]][which(rownames(core_haplotype_abun[[phylotype]]) %in% intersecting_samples), ]
  
  sample_concordance[gene, "RVadj"] <-   MatrixCorrelation::RVadj(as.matrix(otu_tab_subset),
                                                                  as.matrix(strain_tab_subset))
  
}

