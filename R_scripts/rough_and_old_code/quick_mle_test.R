
test_genes_and_num_strains <- read.table(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/inferred_num_strains_pandora_w_struct.tsv",
                                         header = FALSE, sep = "\t", stringsAsFactors = FALSE)

test_genes_and_num_strains <- test_genes_and_num_strains[-which(test_genes_and_num_strains$V2 == 1), ]
tmp <- list()

for (i in 1:nrow(test_genes_and_num_strains)) {
  tmp[[i]] <- test(gene = test_genes_and_num_strains$V1[i], num_orig_haplotypes = test_genes_and_num_strains$V2[i])
}

tmp <- unlist(tmp)

test <- function(gene, num_orig_haplotypes) {

  gene_samples_path <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/prepped_input/", 
                             gene, "_samples.txt", sep = "")
  gene_samples <- read.table(file = gene_samples_path, stringsAsFactors = FALSE)$V1
  
  phylotype <- strsplit(gene, split = "_")[[1]][1]
  phylotype_strains <- core_haplotype_abun[[phylotype]]
  
  otu_tab_path <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/output/pandora_prepped_w_struct/otu_tables/otu_table.",
                        gene, ".", as.character(num_orig_haplotypes), ".txt", sep = "")
  otu_tab <- read.table(file = otu_tab_path, header = FALSE, skip = 1, sep = "\t")
  otu_tab[otu_tab < 0.01] <- 0
  
  rownames(otu_tab) <- gene_samples
  
  # Remove haplotypes that never intersect with any strains.
  gene_and_strain_samples <- gene_samples[which(gene_samples %in% rownames(phylotype_strains))]
  intersecting_haplotypes_i <- which(colSums(otu_tab[gene_and_strain_samples, ]) > 0)
  otu_tab <- otu_tab[, intersecting_haplotypes_i, drop = FALSE]
  
  num_haplotypes <- ncol(otu_tab)
  
  # Add missing samples to the strain and OTU tables (with only abundances given by the pseudocounts below).
  all_samples_to_test <- unique(c(gene_samples, rownames(phylotype_strains)))
  
  gene_samples_with_no_strain <- gene_samples[which(! gene_samples %in% rownames(phylotype_strains))]
  if (length(gene_samples_with_no_strain) > 0) {
    phylotype_strains[gene_samples_with_no_strain, ] <- 0
  }
  
  strain_samples_with_no_gene <- rownames(phylotype_strains)[which(! rownames(phylotype_strains) %in% gene_samples)]
  if (length(strain_samples_with_no_gene) > 0) {
    otu_tab[strain_samples_with_no_gene, ] <- 0
  }
  
  # For each strain, figure out which haplotypes it could ostensibly encode
  # The clearest cases to exclude are those where the strain and haplotype
  # never co-occur in the same sample. Can also ignore strain / haplotype
  # possibilities that don't co-occur at least 75% of the time when the strain is there.
  
  possible_strain_haplotypes <- list()
  
  strains2ignore <- c()
  
  for (strain_i in 1:ncol(phylotype_strains)) {
    
    haplotypes_to_add <- c(1:ncol(otu_tab), ncol(otu_tab) + 1)
    
    samples_w_strain <- rownames(phylotype_strains)[which(phylotype_strains[, strain_i] > 0)]
    
    sample_haplotype_to_ignore <- which(colSums(otu_tab[samples_w_strain, , drop = FALSE] > 0) / length(samples_w_strain) < 0.5)
    
    if (length(sample_haplotype_to_ignore) > 0) {
      if (length(haplotypes_to_add[-sample_haplotype_to_ignore]) == 1 && haplotypes_to_add[-sample_haplotype_to_ignore] == ncol(otu_tab) + 1) {
        strains2ignore <- c(strains2ignore, strain_i)
      }
      
      possible_strain_haplotypes[[strain_i]] <- haplotypes_to_add[-sample_haplotype_to_ignore]
    } else {
      possible_strain_haplotypes[[strain_i]] <- haplotypes_to_add
    }
  }
  
  if (length(strains2ignore) > 0) {
    possible_strain_haplotypes <- possible_strain_haplotypes[-strains2ignore]
    phylotype_strains <- phylotype_strains[, -strains2ignore]
  }

  if (length(possible_strain_haplotypes) > 0) {
    return(prod(sapply(possible_strain_haplotypes, length)))
  } else {
    return(NA)
  }
   
}






max_strain_haplotype_intersect_prop <- function(gene, num_orig_haplotypes) {
  
  gene_samples_path <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/prepped_input/", 
                             gene, "_samples.txt", sep = "")
  gene_samples <- read.table(file = gene_samples_path, stringsAsFactors = FALSE)$V1
  
  phylotype <- strsplit(gene, split = "_")[[1]][1]
  phylotype_strains <- core_haplotype_abun[[phylotype]]
  
  otu_tab_path <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/output/pandora_prepped_w_struct/otu_tables/otu_table.",
                        gene, ".", as.character(num_orig_haplotypes), ".txt", sep = "")
  otu_tab <- read.table(file = otu_tab_path, header = FALSE, skip = 1, sep = "\t")
  otu_tab[otu_tab < 0.01] <- 0
  
  rownames(otu_tab) <- gene_samples
  
  # Remove haplotypes that never intersect with any strains.
  gene_and_strain_samples <- gene_samples[which(gene_samples %in% rownames(phylotype_strains))]
  intersecting_haplotypes_i <- which(colSums(otu_tab[gene_and_strain_samples, ]) > 0)
  otu_tab <- otu_tab[, intersecting_haplotypes_i, drop = FALSE]
  
  num_haplotypes <- ncol(otu_tab)
  
  # Add missing samples to the strain and OTU tables (with only abundances given by the pseudocounts below).
  all_samples_to_test <- unique(c(gene_samples, rownames(phylotype_strains)))
  
  gene_samples_with_no_strain <- gene_samples[which(! gene_samples %in% rownames(phylotype_strains))]
  if (length(gene_samples_with_no_strain) > 0) {
    phylotype_strains[gene_samples_with_no_strain, ] <- 0
  }
  
  strain_samples_with_no_gene <- rownames(phylotype_strains)[which(! rownames(phylotype_strains) %in% gene_samples)]
  if (length(strain_samples_with_no_gene) > 0) {
    otu_tab[strain_samples_with_no_gene, ] <- 0
  }
  
  # For each strain, figure out which haplotypes it could ostensibly encode
  # The clearest cases to exclude are those where the strain and haplotype
  # never co-occur in the same sample. Can also ignore strain / haplotype
  # possibilities that don't co-occur at least 75% of the time when the strain is there.
  
  possible_strain_haplotypes <- list()
  
  prop_intersect <- c()
  
  for (strain_i in 1:ncol(phylotype_strains)) {
    
    haplotypes_to_add <- c(1:ncol(otu_tab), ncol(otu_tab) + 1)
    
    samples_w_strain <- rownames(phylotype_strains)[which(phylotype_strains[, strain_i] > 0)]
    
    if (length(samples_w_strain) == 0) { print(gene) }
    
    prop_intersect <- c(prop_intersect, colSums(otu_tab[samples_w_strain, , drop = FALSE] > 0) / length(samples_w_strain))
    
    
  }
  
  return(max(prop_intersect))
  
}


test_genes_and_num_strains <- read.table(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/inferred_num_strains_pandora_w_struct.tsv",
                                         header = FALSE, sep = "\t", stringsAsFactors = FALSE)

test_genes_and_num_strains <- test_genes_and_num_strains[-which(test_genes_and_num_strains$V2 == 1), ]
tmp <- list()

for (i in 1:nrow(test_genes_and_num_strains)) {
  tmp[[i]] <- max_strain_haplotype_intersect_prop(gene = test_genes_and_num_strains$V1[i], num_orig_haplotypes = test_genes_and_num_strains$V2[i])
}

tmp <- unlist(tmp)

max_intersect_df <- data.frame(max_intersect = tmp, gene = test_genes_and_num_strains$V1)

max_intersect_df$phylotype <- sapply(strsplit(max_intersect_df$gene, "_"), function(x) { x[1] })

library(ggplot2)

ggplot(data = max_intersect_df, aes(x = max_intersect, fill = phylotype)) +
       geom_histogram() +
       facet_grid(. ~ phylotype)

