rm(list = ls(all.names = TRUE))

setwd("honey_bee_pangenome/data/Ellegaard_Carrie_pangenome_mapping/")

Snodgrassella_subfam80_info <- read.table("subfam_info/Snodgrassella/snod_subfamilies_list_80.txt",
                                          stringsAsFactors = FALSE, sep = "\t",
                                          col.names = paste0("V", seq_len(57 + 4)), fill = TRUE)

determine_unique_genomes <- function(info_tab) {

  unique_genomes <- c()

  for(row_i in 1:nrow(info_tab)) {
      
    tmp_unique <- sapply(info_tab[row_i, 5:ncol(info_tab)], substr, 1, 4)  
    names(tmp_unique) <- NULL
    
    unique_genomes <- c(unique_genomes, tmp_unique)
      
    if(length(which(duplicated(unique_genomes))) > 0) {
      unique_genomes <- unique_genomes[-which(duplicated(unique_genomes))]
    }
  }

  if ("" %in% unique_genomes) {
    unique_genomes <- unique_genomes[-which(unique_genomes == "")] 
  }
  
  return(unique_genomes)
}

Snodgrassella_subfam80_unique_genomes <- determine_unique_genomes(info_tab = Snodgrassella_subfam80_info)

Snodgrassella_subfam80_info_core <- Snodgrassella_subfam80_info[which(Snodgrassella_subfam80_info$V3 == 57), ]

Snodgrassella_subfam80_info_core$core <- FALSE

for(row_i in 1:nrow(Snodgrassella_subfam80_info_core)) {
  Snodgrassella_subfam80_info_core[row_i, 5:ncol(Snodgrassella_subfam80_info_core)] <- sapply(Snodgrassella_subfam80_info_core[row_i, 5:ncol(Snodgrassella_subfam80_info_core)], substr, 1, 4)
  
  if(length(which(Snodgrassella_subfam80_unique_genomes %in% Snodgrassella_subfam80_info_core[row_i, 5:ncol(Snodgrassella_subfam80_info_core)])) == 57) {
    Snodgrassella_subfam80_info_core[row_i, "core"] <- TRUE
  }
}

Snodgrassella_subfam80_info_core <- Snodgrassella_subfam80_info_core[which(Snodgrassella_subfam80_info_core$core == TRUE), ]
Snodgrassella_subfam80_info_core_geneids <- Snodgrassella_subfam80_info_core$V1
Snodgrassella_subfam80_info_core_geneids <- paste("Snod", Snodgrassella_subfam80_info_core_geneids, "c1", sep = "_")


trimmed_SNVs <- read.table("instrain_out/subfam80/subfam80_Snodgrassella_all_unique_SNVs.trimmedgenes.bed.gz",
                           header = FALSE, sep = "\t", stringsAsFactors = FALSE)
trimmed_SNVs <- paste(trimmed_SNVs$V1, trimmed_SNVs$V2, sep = "|")

instrain_SNV <- read.table("instrain_out/subfam80/Snodgrassella/SRR7287205_SNVs.tsv",
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(instrain_SNV) <- paste(instrain_SNV$scaffold, instrain_SNV$position, sep = "|")


core_instrain_SNV <- instrain_SNV[which(instrain_SNV$scaffold %in% Snodgrassella_subfam80_info_core_geneids), ]
noncore_instrain_SNV <- instrain_SNV[which(! instrain_SNV$scaffold %in% Snodgrassella_subfam80_info_core_geneids), ]

par(mfrow = c(2, 1))
core_SNV_depth <- c(core_instrain_SNV$A, core_instrain_SNV$C, core_instrain_SNV$T, core_instrain_SNV$G)
core_SNV_depth <- core_SNV_depth[-which(core_SNV_depth == 0)]
hist(core_SNV_depth, breaks = 1000, xlim=c(0, 100))

noncore_SNV_depth <- c(noncore_instrain_SNV$A, noncore_instrain_SNV$C, noncore_instrain_SNV$T, noncore_instrain_SNV$G)
noncore_SNV_depth <- noncore_SNV_depth[-which(noncore_SNV_depth == 0)]
hist(noncore_SNV_depth, breaks = 1000, xlim=c(0, 100))

instrain_SNV_trimmed <- instrain_SNV[which(rownames(instrain_SNV) %in% trimmed_SNVs), ]

all_depth_trimmed <- c(instrain_SNV_trimmed$A, instrain_SNV_trimmed$C,
                       instrain_SNV_trimmed$T, instrain_SNV_trimmed$G)
all_depth_trimmed <- all_depth_trimmed[-which(all_depth_trimmed == 0)]

hist(all_depth_trimmed, breaks = 1000)



par(mfrow = c(2, 1))
all_depth <- c(instrain_SNV$A, instrain_SNV$C, instrain_SNV$T, instrain_SNV$G)
all_depth <- all_depth[-which(all_depth == 0)]

hist(all_depth, breaks = 1000)


instrain_SNV_trimmed <- instrain_SNV[which(rownames(instrain_SNV) %in% trimmed_SNVs), ]

all_depth_trimmed <- c(instrain_SNV_trimmed$A, instrain_SNV_trimmed$C,
                       instrain_SNV_trimmed$T, instrain_SNV_trimmed$G)
all_depth_trimmed <- all_depth_trimmed[-which(all_depth_trimmed == 0)]

hist(all_depth_trimmed, breaks = 1000)




gene_depth_untrimmed <- read.table("gene_depth_coverage/subfam80_gene/SRR7287194.Snodgrassella.subfam80.gene.coverage.bed.gz",
                                   header = FALSE, sep = "\t", stringsAsFactors = FALSE)

gene_depth_trimmed <- read.table("gene_depth_coverage/subfam80_gene.trimmed/SRR7287194.Snodgrassella.subfam80.gene.trimmed.coverage.bed.gz",
                                 header = FALSE, sep = "\t", stringsAsFactors = FALSE)

colnames(gene_depth_untrimmed) <- c("scaffold", "start", "stop", "name", "num_reads", "num_covered_bases", "gene_length", "fraction_covered_based")
colnames(gene_depth_trimmed) <- c("scaffold", "start", "stop", "name", "num_reads", "num_covered_bases", "gene_length", "fraction_covered_based")

duplicated_genes <- gene_depth_untrimmed$scaffold[which(duplicated(gene_depth_untrimmed$scaffold))]
gene_depth_untrimmed <- gene_depth_untrimmed[-which(gene_depth_untrimmed$scaffold %in% duplicated_genes), ]
gene_depth_trimmed <- gene_depth_trimmed[-which(gene_depth_trimmed$scaffold %in% duplicated_genes), ]

rownames(gene_depth_untrimmed) <- gene_depth_untrimmed$scaffold
rownames(gene_depth_trimmed) <- gene_depth_trimmed$scaffold

gene_depth_untrimmed <- gene_depth_untrimmed[-which(gene_depth_untrimmed$num_reads == 0), ]

intersecting_genes <- rownames(gene_depth_untrimmed)[which(rownames(gene_depth_untrimmed) %in% rownames(gene_depth_trimmed))]

gene_depth_untrimmed <- gene_depth_untrimmed[intersecting_genes, ]
gene_depth_trimmed <- gene_depth_trimmed[intersecting_genes, ]

gene_depth_untrimmed$mean_depth <- gene_depth_untrimmed$num_reads / gene_depth_untrimmed$gene_length
gene_depth_trimmed$mean_depth <- gene_depth_trimmed$num_reads / gene_depth_trimmed$gene_length

hist(gene_depth_untrimmed$mean_depth, breaks = 1000, xlim = c(0, 0.6))
hist(gene_depth_trimmed$mean_depth, breaks = 1000, xlim = c(0, 0.6))

plot(gene_depth_untrimmed$mean_depth, gene_depth_trimmed$mean_depth)
