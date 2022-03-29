### Parse FastANI of final genomes and quickly check that they make sense visually.

rm(list = ls(all.names = TRUE))

library(ape)
library(ggplot2)
library(ggbeeswarm)
library(ggtree)

setwd("/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes_draft_ANI")

generate_fastani_dist_df <- function(in_df) {
  
  all_g <- unique(c(in_df$G1, in_df$G2))
  
  out_df <- data.frame(matrix(NA, nrow = length(all_g), ncol = length(all_g)))
  
  rownames(out_df) <- all_g
  colnames(out_df) <- all_g
  
  for (g1 in all_g) {
    
    g1_i <- which(all_g == g1)
    
    if (g1_i == length(all_g)) { next }
    
    remaining_g <- all_g[(g1_i + 1):length(all_g)]
    
    for (g2 in remaining_g) {
      
      matching_rows <- c(which(in_df$G1 == g1 & in_df$G2 == g2), which(in_df$G1 == g2 & in_df$G2 == g1))
      
      if (length(matching_rows) == 0) { next }
      
      ANI <- in_df[matching_rows, "ani"][1]
      out_df[g1, g2] <- ANI
      out_df[g2, g1] <- ANI
      
    }
    
  }
  
  return(out_df)
  
}

fastani_out <- read.table("final_genomes_fastani/fastani_output.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

fastani_out$V1 <- gsub(".fna", "", fastani_out$V1)
fastani_out$V2 <- gsub(".fna", "", fastani_out$V2)
fastani_out$V1 <- gsub("../../highqual_genomes/", "", fastani_out$V1)
fastani_out$V2 <- gsub("../../highqual_genomes/", "", fastani_out$V2)

colnames(fastani_out) <- c("file1", "file2", "ani", "num_aligned", "num_fragments")

fastani_out$sp1 <- gsub("/.*$", "", fastani_out$file1)
fastani_out$id1 <- gsub("^.*/", "", fastani_out$file1)

fastani_out$sp2 <- gsub("/.*$", "", fastani_out$file2)
fastani_out$id2 <- gsub("^.*/", "", fastani_out$file2)


species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species.txt", stringsAsFactors = FALSE)$V1

fastani_by_species <- list()

for (s in species) {
  
  fastani_file1_match_i <- grep(s, fastani_out$file1)
  fastani_file2_match_i <- grep(s, fastani_out$file2)
  
  fastani_file_matches <- sort(unique(c(fastani_file1_match_i, fastani_file2_match_i)))
  
  if (nrow(fastani_out[fastani_file_matches,]) == 0) { next }
  
  fastani_by_species[[s]] <- fastani_out[fastani_file_matches,]
  
  fastani_by_species[[s]]$sp1 <- gsub("_", " ", fastani_by_species[[s]]$sp1)
  fastani_by_species[[s]]$sp2 <- gsub("_", " ", fastani_by_species[[s]]$sp2)
  
  fastani_by_species[[s]]$G1 <- paste(fastani_by_species[[s]]$sp1, fastani_by_species[[s]]$id1, sep = "|")
  fastani_by_species[[s]]$G2 <- paste(fastani_by_species[[s]]$sp2, fastani_by_species[[s]]$id2, sep = "|")
  
  fastani_by_species[[s]]$G_combined <- NA
  fastani_by_species[[s]]$sp_combined <- NA
  
  for (i in 1:nrow(fastani_by_species[[s]])) {
   G_combined <- sort(c(fastani_by_species[[s]][i, "G1"], fastani_by_species[[s]][i, "G2"]))
   fastani_by_species[[s]]$G_combined[i] <- paste(G_combined, collapse = " ")
   
   sp_combined <- sort(c(fastani_by_species[[s]][i, "sp1"], fastani_by_species[[s]][i, "sp2"]))
   fastani_by_species[[s]]$sp_combined[i] <- paste(sp_combined, collapse = " vs. ")
   
  }
  
  print(c(s, length(unique(c(fastani_by_species[[s]]$id1, fastani_by_species[[s]]$id2)))))
  
}




### Apilactobacillus_apinorum
Apilactobacillus_apinorum <- fastani_by_species$Apilactobacillus_apinorum
Apilactobacillus_apinorum <- Apilactobacillus_apinorum[-which(duplicated(Apilactobacillus_apinorum$G_combined)), ]
Apilactobacillus_apinorum <- Apilactobacillus_apinorum[-which(Apilactobacillus_apinorum$id1 == Apilactobacillus_apinorum$id2), ]
ggplot(data = Apilactobacillus_apinorum, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

### Apilactobacillus_kunkeei
Apilactobacillus_kunkeei <- fastani_by_species$Apilactobacillus_kunkeei
Apilactobacillus_kunkeei <- Apilactobacillus_kunkeei[-which(duplicated(Apilactobacillus_kunkeei$G_combined)), ]
Apilactobacillus_kunkeei <- Apilactobacillus_kunkeei[-which(Apilactobacillus_kunkeei$id1 == Apilactobacillus_kunkeei$id2), ]
ggplot(data = Apilactobacillus_kunkeei, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

### Bartonella_apis
Bartonella_apis <- fastani_by_species$Bartonella_apis
Bartonella_apis <- Bartonella_apis[-which(duplicated(Bartonella_apis$G_combined)), ]
Bartonella_apis <- Bartonella_apis[-which(Bartonella_apis$id1 == Bartonella_apis$id2), ]
ggplot(data = Bartonella_apis, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


### Bifidobacterium_asteroides

Bifidobacterium_asteroides <- fastani_by_species$Bifidobacterium_asteroides
Bifidobacterium_asteroides <- Bifidobacterium_asteroides[-which(duplicated(Bifidobacterium_asteroides$G_combined)), ]
Bifidobacterium_asteroides <- Bifidobacterium_asteroides[-which(Bifidobacterium_asteroides$id1 == Bifidobacterium_asteroides$id2), ]
ggplot(data = Bifidobacterium_asteroides, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


### Bifidobacterium_coryneforme_indicum
Bifidobacterium_coryneforme_indicum <- fastani_by_species$Bifidobacterium_coryneforme_indicum
Bifidobacterium_coryneforme_indicum <- Bifidobacterium_coryneforme_indicum[-which(duplicated(Bifidobacterium_coryneforme_indicum$G_combined)), ]
Bifidobacterium_coryneforme_indicum <- Bifidobacterium_coryneforme_indicum[-which(Bifidobacterium_coryneforme_indicum$id1 == Bifidobacterium_coryneforme_indicum$id2), ]
ggplot(data = Bifidobacterium_coryneforme_indicum, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


### Bombella_apis
Bombella_apis <- fastani_by_species$Bombella_apis
Bombella_apis <- Bombella_apis[-which(duplicated(Bombella_apis$G_combined)), ]
Bombella_apis <- Bombella_apis[-which(Bombella_apis$id1 == Bombella_apis$id2), ]
ggplot(data = Bombella_apis, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


### Bombella_sp
Bombella_sp <- fastani_by_species$Bombella_sp
Bombella_sp <- Bombella_sp[-which(duplicated(Bombella_sp$G_combined)), ]
Bombella_sp <- Bombella_sp[-which(Bombella_sp$id1 == Bombella_sp$id2), ]
ggplot(data = Bombella_sp, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")



### Bombilactobacillus_mellifer
Bombilactobacillus_mellifer <- fastani_by_species$Bombilactobacillus_mellifer
Bombilactobacillus_mellifer <- Bombilactobacillus_mellifer[-which(duplicated(Bombilactobacillus_mellifer$G_combined)), ]
Bombilactobacillus_mellifer <- Bombilactobacillus_mellifer[-which(Bombilactobacillus_mellifer$id1 == Bombilactobacillus_mellifer$id2), ]
ggplot(data = Bombilactobacillus_mellifer, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


### Bombilactobacillus_mellis
Bombilactobacillus_mellis <- fastani_by_species$Bombilactobacillus_mellis
Bombilactobacillus_mellis <- Bombilactobacillus_mellis[-which(duplicated(Bombilactobacillus_mellis$G_combined)), ]
Bombilactobacillus_mellis <- Bombilactobacillus_mellis[-which(Bombilactobacillus_mellis$id1 == Bombilactobacillus_mellis$id2), ]
ggplot(data = Bombilactobacillus_mellis, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


### Commensalibacter_sp
Commensalibacter_sp <- fastani_by_species$Commensalibacter_sp
Commensalibacter_sp <- Commensalibacter_sp[-which(duplicated(Commensalibacter_sp$G_combined)), ]
Commensalibacter_sp <- Commensalibacter_sp[-which(Commensalibacter_sp$id1 == Commensalibacter_sp$id2), ]
ggplot(data = Commensalibacter_sp, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


### Frischella_perrara
Frischella_perrara <- fastani_by_species$Frischella_perrara
Frischella_perrara <- Frischella_perrara[-which(duplicated(Frischella_perrara$G_combined)), ]
Frischella_perrara <- Frischella_perrara[-which(Frischella_perrara$id1 == Frischella_perrara$id2), ]
ggplot(data = Frischella_perrara, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


### Gilliamella_apicola
Gilliamella_apicola <- fastani_by_species$Gilliamella_apicola
Gilliamella_apicola <- Gilliamella_apicola[-which(duplicated(Gilliamella_apicola$G_combined)), ]
Gilliamella_apicola <- Gilliamella_apicola[-which(Gilliamella_apicola$id1 == Gilliamella_apicola$id2), ]
ggplot(data = Gilliamella_apicola, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


### Gilliamella_apis
Gilliamella_apis <- fastani_by_species$Gilliamella_apis
Gilliamella_apis <- Gilliamella_apis[-which(duplicated(Gilliamella_apis$G_combined)), ]
Gilliamella_apis <- Gilliamella_apis[-which(Gilliamella_apis$id1 == Gilliamella_apis$id2), ]
ggplot(data = Gilliamella_apis, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


### Gilliamella_sp
Gilliamella_sp <- fastani_by_species$Gilliamella_sp
Gilliamella_sp <- Gilliamella_sp[-which(duplicated(Gilliamella_sp$G_combined)), ]
Gilliamella_sp <- Gilliamella_sp[-which(Gilliamella_sp$id1 == Gilliamella_sp$id2), ]
ggplot(data = Gilliamella_sp, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


### Lactobacillus_apis
Lactobacillus_apis <- fastani_by_species$Lactobacillus_apis
Lactobacillus_apis <- Lactobacillus_apis[-which(duplicated(Lactobacillus_apis$G_combined)), ]
Lactobacillus_apis <- Lactobacillus_apis[-which(Lactobacillus_apis$id1 == Lactobacillus_apis$id2), ]
ggplot(data = Lactobacillus_apis, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


### Lactobacillus_helsingborgensis
Lactobacillus_helsingborgensis <- fastani_by_species$Lactobacillus_helsingborgensis
Lactobacillus_helsingborgensis <- Lactobacillus_helsingborgensis[-which(duplicated(Lactobacillus_helsingborgensis$G_combined)), ]
Lactobacillus_helsingborgensis <- Lactobacillus_helsingborgensis[-which(Lactobacillus_helsingborgensis$id1 == Lactobacillus_helsingborgensis$id2), ]
ggplot(data = Lactobacillus_helsingborgensis, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

### Lactobacillus_kimbladii
Lactobacillus_kimbladii <- fastani_by_species$Lactobacillus_kimbladii
Lactobacillus_kimbladii <- Lactobacillus_kimbladii[-which(duplicated(Lactobacillus_kimbladii$G_combined)), ]
Lactobacillus_kimbladii <- Lactobacillus_kimbladii[-which(Lactobacillus_kimbladii$id1 == Lactobacillus_kimbladii$id2), ]
ggplot(data = Lactobacillus_kimbladii, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

### Lactobacillus_kullabergensis
Lactobacillus_kullabergensis <- fastani_by_species$Lactobacillus_kullabergensis
Lactobacillus_kullabergensis <- Lactobacillus_kullabergensis[-which(duplicated(Lactobacillus_kullabergensis$G_combined)), ]
Lactobacillus_kullabergensis <- Lactobacillus_kullabergensis[-which(Lactobacillus_kullabergensis$id1 == Lactobacillus_kullabergensis$id2), ]
ggplot(data = Lactobacillus_kullabergensis, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

### Lactobacillus_melliventris
Lactobacillus_melliventris <- fastani_by_species$Lactobacillus_melliventris
Lactobacillus_melliventris <- Lactobacillus_melliventris[-which(duplicated(Lactobacillus_melliventris$G_combined)), ]
Lactobacillus_melliventris <- Lactobacillus_melliventris[-which(Lactobacillus_melliventris$id1 == Lactobacillus_melliventris$id2), ]
ggplot(data = Lactobacillus_melliventris, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

### Serratia_marcescens
Serratia_marcescens <- fastani_by_species$Serratia_marcescens
Serratia_marcescens <- Serratia_marcescens[-which(duplicated(Serratia_marcescens$G_combined)), ]
Serratia_marcescens <- Serratia_marcescens[-which(Serratia_marcescens$id1 == Serratia_marcescens$id2), ]
ggplot(data = Serratia_marcescens, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

### Snodgrassella_alvi
Snodgrassella_alvi <- fastani_by_species$Snodgrassella_alvi
Snodgrassella_alvi <- Snodgrassella_alvi[-which(duplicated(Snodgrassella_alvi$G_combined)), ]
Snodgrassella_alvi <- Snodgrassella_alvi[-which(Snodgrassella_alvi$id1 == Snodgrassella_alvi$id2), ]
ggplot(data = Snodgrassella_alvi, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")




