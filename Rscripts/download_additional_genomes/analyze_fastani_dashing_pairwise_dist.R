### Compare FastANI / Dashing distances between all pairwise high quality genomes.

rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
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

fastani_out <- read.table("fastani/fastani_output.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

fastani_out$V1 <- gsub(".fna", "", fastani_out$V1)
fastani_out$V2 <- gsub(".fna", "", fastani_out$V2)
fastani_out$V1 <- gsub("../../highqual_genomes_draft/", "", fastani_out$V1)
fastani_out$V2 <- gsub("../../highqual_genomes_draft/", "", fastani_out$V2)

colnames(fastani_out) <- c("file1", "file2", "ani", "num_aligned", "num_fragments")

fastani_out$sp1 <- gsub("/.*$", "", fastani_out$file1)
fastani_out$id1 <- gsub("^.*/", "", fastani_out$file1)

fastani_out$sp2 <- gsub("/.*$", "", fastani_out$file2)
fastani_out$id2 <- gsub("^.*/", "", fastani_out$file2)


species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/raw_species.txt", stringsAsFactors = FALSE)$V1

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


### Apilactobacillus apinoruim vs kunkeei

Apilactobacillus_kunkeei_fastani_dist_df <- generate_fastani_dist_df(fastani_by_species$Apilactobacillus_kunkeei)

Apilactobacillus_kunkeei_fastani_phylo <- ape::njs(as.dist(100 - Apilactobacillus_kunkeei_fastani_dist_df))

Apilactobacillus_kunkeei_final_set <- read.table("../highqual_genomes_draft/Apilactobacillus_kunkeei.txt", stringsAsFactors = FALSE)$V1
Apilactobacillus_apinorum_final_set <- read.table("../highqual_genomes_draft/Apilactobacillus_apinorum.txt", stringsAsFactors = FALSE)$V1

Apilactobacillus_kunkeei_final_set <- c(Apilactobacillus_kunkeei_final_set, "GCF_016101275.1")

write.table(x = Apilactobacillus_kunkeei_final_set,
            file = "../highqual_genomes/Apilactobacillus_kunkeei.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(x = Apilactobacillus_apinorum_final_set,
            file = "../highqual_genomes/Apilactobacillus_apinorum.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

### Bartonella_apis

Bartonella_apis_out <- fastani_by_species$Bartonella_apis
Bartonella_apis_out <- Bartonella_apis_out[-which(duplicated(Bartonella_apis_out$G_combined)), ]
Bartonella_apis_out <- Bartonella_apis_out[-which(Bartonella_apis_out$id1 == Bartonella_apis_out$id2), ]
ggplot(data = Bartonella_apis_out, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

Bartonella_apis_fastani_dist_df <- generate_fastani_dist_df(fastani_by_species$Bartonella_apis)
Bartonella_apis_fastani_phylo <- ape::njs(as.dist(100 - Bartonella_apis_fastani_dist_df))
plot(Bartonella_apis_fastani_phylo)

Bartonella_apis_final_set <- read.table("../highqual_genomes_draft/Bartonella_apis.txt", stringsAsFactors = FALSE)$V1

write.table(x = Bartonella_apis_final_set,
            file = "../highqual_genomes/Bartonella_apis.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

### All Bifidobacterium

all_bifido_out <- rbind(fastani_by_species$Bifidobacterium_asteroides, fastani_by_species$Bifidobacterium_coryneforme)
all_bifido_out <- rbind(all_bifido_out, fastani_by_species$Bifidobacterium_indicum)
all_bifido_out <- all_bifido_out[-which(duplicated(all_bifido_out$G_combined)), ]
all_bifido_out <- all_bifido_out[-which(all_bifido_out$id1 == all_bifido_out$id2), ]
ggplot(data = all_bifido_out, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


all_bifido_dist_df <- generate_fastani_dist_df(all_bifido_out)
all_bifido_fastani_phylo <- ape::njs(as.dist(100 - all_bifido_dist_df))
plot(all_bifido_fastani_phylo)


Bifidobacterium_asteroides_final_set <- read.table("../highqual_genomes_draft/Bifidobacterium_asteroides.txt", stringsAsFactors = FALSE)$V1
write.table(x = Bifidobacterium_asteroides_final_set,
            file = "../highqual_genomes/Bifidobacterium_asteroides.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

Bifidobacterium_coryneforme_orig_set <- read.table("../highqual_genomes_draft/Bifidobacterium_coryneforme.txt", stringsAsFactors = FALSE)$V1
Bifidobacterium_indicum_orig_set <- read.table("../highqual_genomes_draft/Bifidobacterium_indicum.txt", stringsAsFactors = FALSE)$V1

Bifidobacterium_coryneforme_indicum_final_set <- c(Bifidobacterium_coryneforme_orig_set, Bifidobacterium_indicum_orig_set)

write.table(x = Bifidobacterium_coryneforme_indicum_final_set,
            file = "../highqual_genomes/Bifidobacterium_coryneforme_indicum.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)



### Bombella apis / Bombella sp / Parasaccharibacter_apium / Saccharibacter_sp

all_bombella_types <- rbind(fastani_by_species$Bombella_apis, fastani_by_species$Bombella_sp)
all_bombella_types <- rbind(all_bombella_types, fastani_by_species$Parasaccharibacter_apium)
all_bombella_types <- rbind(all_bombella_types, fastani_by_species$Saccharibacter_sp)
all_bombella_types <- all_bombella_types[-which(duplicated(all_bombella_types$G_combined)), ]
all_bombella_types <- all_bombella_types[-which(all_bombella_types$id1 == all_bombella_types$id2), ]
ggplot(data = all_bombella_types, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")



all_bombella_types_dist_df <- generate_fastani_dist_df(all_bombella_types)

all_bombella_types_fastani_phylo <- ape::njs(as.dist(100 - all_bombella_types_dist_df))

plot(all_bombella_types_fastani_phylo)


Bombella_sp_final_set <- c("GCF_002592045.1", "GCF_009725755.1", "GCF_009725845.1")

Bombella_apis_orig_set <- read.table("../highqual_genomes_draft/Bombella_apis.txt", stringsAsFactors = FALSE)$V1
Bombella_sp_orig_set <- read.table("../highqual_genomes_draft/Bombella_sp.txt", stringsAsFactors = FALSE)$V1
Parasaccharibacter_apium_orig_set <- read.table("../highqual_genomes_draft/Parasaccharibacter_apium.txt", stringsAsFactors = FALSE)$V1
Saccharibacter_sp_orig_set <- read.table("../highqual_genomes_draft/Saccharibacter_sp.txt", stringsAsFactors = FALSE)$V1

Bombella_apis_final_set <- c(Bombella_apis_orig_set, Bombella_sp_orig_set, Parasaccharibacter_apium_orig_set, Saccharibacter_sp_orig_set)
Bombella_apis_final_set <- Bombella_apis_final_set[-which(Bombella_apis_final_set %in% Bombella_sp_final_set)]

write.table(x = Bombella_sp_final_set,
            file = "../highqual_genomes/Bombella_sp.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(x = Bombella_apis_final_set,
            file = "../highqual_genomes/Bombella_apis.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)


### Bombilactobacillus mellis / mellifer
bombilactobacillus <- rbind(fastani_by_species$Bombilactobacillus_mellifer, fastani_by_species$Bombilactobacillus_mellis)
bombilactobacillus <- bombilactobacillus[-which(duplicated(bombilactobacillus$G_combined)), ]
bombilactobacillus <- bombilactobacillus[-which(bombilactobacillus$id1 == bombilactobacillus$id2), ]
ggplot(data = bombilactobacillus, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

bombilactobacillus_dist_df <- generate_fastani_dist_df(bombilactobacillus)
bombilactobacillus_fastani_phylo <- ape::njs(as.dist(100 - bombilactobacillus_dist_df))
plot(bombilactobacillus_fastani_phylo)


Bombilactobacillus_mellifer_orig_set <- read.table("../highqual_genomes_draft/Bombilactobacillus_mellifer.txt", stringsAsFactors = FALSE)$V1
Bombilactobacillus_mellis_orig_set <- read.table("../highqual_genomes_draft/Bombilactobacillus_mellis.txt", stringsAsFactors = FALSE)$V1

Bombilactobacillus_mellis_final_set <- c(Bombilactobacillus_mellis_orig_set,
                                           "GCF_016100965.1", "GCF_016100985.1", "GCF_016101055.1",
                                           "GCF_016101025.1", "GCF_016101075.1", "GCF_016101045.1")

Bombilactobacillus_mellifer_final_set <- c(Bombilactobacillus_mellifer_orig_set, "GCF_016102125.1")

write.table(x = Bombilactobacillus_mellifer_final_set,
            file = "../highqual_genomes/Bombilactobacillus_mellifer.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(x = Bombilactobacillus_mellis_final_set,
            file = "../highqual_genomes/Bombilactobacillus_mellis.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)



### Commensalibacter_sp
Commensalibacter <- fastani_by_species$Commensalibacter_sp
Commensalibacter <- Commensalibacter[-which(duplicated(Commensalibacter$G_combined)), ]
Commensalibacter <- Commensalibacter[-which(Commensalibacter$id1 == Commensalibacter$id2), ]
ggplot(data = Commensalibacter, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

Commensalibacter_dist_df <- generate_fastani_dist_df(Commensalibacter)
Commensalibacter_fastani_phylo <- ape::njs(as.dist(100 - Commensalibacter_dist_df))
plot(Commensalibacter_fastani_phylo)

Commensalibacter_sp_final_set <- read.table("../highqual_genomes_draft/Commensalibacter_sp.txt", stringsAsFactors = FALSE)$V1
write.table(x = Commensalibacter_sp_final_set,
            file = "../highqual_genomes/Commensalibacter_sp.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)



#### Frischella_perrara / Gilliamella apis / Gilliamella apicola


gamma_proteo <- rbind(fastani_by_species$Frischella_perrara, fastani_by_species$Gilliamella_apicola)
gamma_proteo <- rbind(gamma_proteo, fastani_by_species$Gilliamella_apis)
gamma_proteo <- gamma_proteo[-which(duplicated(gamma_proteo$G_combined)), ]
gamma_proteo <- gamma_proteo[-which(gamma_proteo$id1 == gamma_proteo$id2), ]
ggplot(data = gamma_proteo, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

Frischella_perrara_final_set <- read.table("../highqual_genomes_draft/Frischella_perrara.txt", stringsAsFactors = FALSE)$V1
write.table(x = Frischella_perrara_final_set,
            file = "../highqual_genomes/Frischella_perrara.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)


gilli <- rbind(fastani_by_species$Gilliamella_apis, fastani_by_species$Gilliamella_apicola)
gilli <- gilli[-which(duplicated(gilli$G_combined)), ]
gilli <- gilli[-which(gilli$id1 == gilli$id2), ]
gilli <- gilli[-which(gilli$sp1 == "Frischella perrara"), ]
ggplot(data = gilli, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

gilli_dist_df <- generate_fastani_dist_df(gilli)

gilli_fastani_phylo <- ape::njs(as.dist(100 - gilli_dist_df))

ggtree(gilli_fastani_phylo) + geom_tiplab()

# Genomes that are major outliers from other Gilliamella apicola
# Gilliamella apicola|GCA_001723875.1
# Gilliamella apicola|GCA_007559165.1
# Gilliamella apicola|GCA_003202915.1
# Gilliamella apicola|GCA_003202655.1

gilli_dist_df["Gilliamella apicola|GCA_001723875.1", grep("Gilliamella apicola", colnames(gilli_dist_df))]
gilli_dist_df["Gilliamella apicola|GCA_001723875.1", grep("Gilliamella apis", colnames(gilli_dist_df))]

gilli_dist_df["Gilliamella apicola|GCA_007559165.1", grep("Gilliamella apicola", colnames(gilli_dist_df))]
gilli_dist_df["Gilliamella apicola|GCA_007559165.1", grep("Gilliamella apis", colnames(gilli_dist_df))]


gilli_dist_df["Gilliamella apicola|GCA_003202915.1", grep("Gilliamella apicola", colnames(gilli_dist_df))]
gilli_dist_df["Gilliamella apicola|GCA_003202915.1", grep("Gilliamella apis", colnames(gilli_dist_df))]

gilli_dist_df["Gilliamella apicola|GCA_003202655.1", grep("Gilliamella apicola", colnames(gilli_dist_df))]
gilli_dist_df["Gilliamella apicola|GCA_003202655.1", grep("Gilliamella apis", colnames(gilli_dist_df))]


gilli_edit <- gilli
gilli_edit$sp1[which(gilli_edit$G1 == "Gilliamella apicola|GCA_001723875.1")] <- "Gilliamella apis"
gilli_edit$sp2[which(gilli_edit$G2 == "Gilliamella apicola|GCA_001723875.1")] <- "Gilliamella apis"

gilli_edit$G1[which(gilli_edit$G1 == "Gilliamella apicola|GCA_001723875.1")] <- "Gilliamella apis|GCA_001723875.1"
gilli_edit$G2[which(gilli_edit$G2 == "Gilliamella apicola|GCA_001723875.1")] <- "Gilliamella apis|GCA_001723875.1"


gilli_edit$sp1[which(gilli_edit$G1 == "Gilliamella apicola|GCA_007559165.1")] <- "Gilliamella apis"
gilli_edit$sp2[which(gilli_edit$G2 == "Gilliamella apicola|GCA_007559165.1")] <- "Gilliamella apis"

gilli_edit$G1[which(gilli_edit$G1 == "Gilliamella apicola|GCA_007559165.1")] <- "Gilliamella apis|GCA_007559165.1"
gilli_edit$G2[which(gilli_edit$G2 == "Gilliamella apicola|GCA_007559165.1")] <- "Gilliamella apis|GCA_007559165.1"


gilli_edit$sp1[which(gilli_edit$G1 == "Gilliamella apicola|GCA_003202915.1")] <- "Gilliamella sp"
gilli_edit$sp2[which(gilli_edit$G2 == "Gilliamella apicola|GCA_003202915.1")] <- "Gilliamella sp"

gilli_edit$G1[which(gilli_edit$G1 == "Gilliamella apicola|GCA_003202915.1")] <- "Gilliamella sp|GCA_003202915.1"
gilli_edit$G2[which(gilli_edit$G2 == "Gilliamella apicola|GCA_003202915.1")] <- "Gilliamella sp|GCA_003202915.1"



gilli_edit$sp1[which(gilli_edit$G1 == "Gilliamella apicola|GCA_003202655.1")] <- "Gilliamella sp"
gilli_edit$sp2[which(gilli_edit$G2 == "Gilliamella apicola|GCA_003202655.1")] <- "Gilliamella sp"

gilli_edit$G1[which(gilli_edit$G1 == "Gilliamella apicola|GCA_003202655.1")] <- "Gilliamella sp|GCA_003202655.1"
gilli_edit$G2[which(gilli_edit$G2 == "Gilliamella apicola|GCA_003202655.1")] <- "Gilliamella sp|GCA_003202655.1"

gilli_edit$sp_combined <- NA

for (i in 1:nrow(gilli_edit)) {
  sp_combined <- sort(c(gilli_edit[i, "sp1"], gilli_edit[i, "sp2"]))
  gilli_edit$sp_combined[i] <- paste(sp_combined, collapse = " vs. ")
}

ggplot(data = gilli_edit, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")




Gilliamella_sp_final_set <- c("GCA_003202915.1", "GCA_003202655.1")

write.table(x = Gilliamella_sp_final_set,
            file = "../highqual_genomes/Gilliamella_sp.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)



Gilliamella_apis_orig_set <- read.table("../highqual_genomes_draft/Gilliamella_apis.txt", stringsAsFactors = FALSE)$V1
Gilliamella_apicola_orig_set <- read.table("../highqual_genomes_draft/Gilliamella_apicola.txt", stringsAsFactors = FALSE)$V1

Gilliamella_apis_final_set <- c(Gilliamella_apis_orig_set, "GCA_007559165.1", "GCA_001723875.1")
write.table(x = Gilliamella_apis_final_set,
            file = "../highqual_genomes/Gilliamella_apis.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

Gilliamella_apicola_final_set <- Gilliamella_apicola_orig_set[-which(Gilliamella_apicola_orig_set %in% c("GCA_007559165.1", "GCA_001723875.1", Gilliamella_sp_final_set))]
write.table(x = Gilliamella_apicola_final_set,
            file = "../highqual_genomes/Gilliamella_apicola.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)


### Lactobacilli

lactobacilli_out <- rbind(fastani_by_species$Lactobacillus_apis, fastani_by_species$Lactobacillus_helsingborgensis)
lactobacilli_out <- rbind(lactobacilli_out, fastani_by_species$Lactobacillus_kimbladii)
lactobacilli_out <- rbind(lactobacilli_out, fastani_by_species$Lactobacillus_kullabergensis)
lactobacilli_out <- rbind(lactobacilli_out, fastani_by_species$Lactobacillus_melliventris)


lactobacilli_out <- lactobacilli_out[-which(duplicated(lactobacilli_out$G_combined)), ]
lactobacilli_out <- lactobacilli_out[-which(lactobacilli_out$id1 == lactobacilli_out$id2), ]


ggplot(data = lactobacilli_out, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")


lactobacilli_dist_df <- generate_fastani_dist_df(lactobacilli_out)
lactobacilli_fastani_phylo <- ape::njs(as.dist(100 - lactobacilli_dist_df))
plot(lactobacilli_fastani_phylo)


Lactobacillus_sp_orig_set <- sort(read.table("../highqual_genomes_draft/Lactobacillus_sp.txt", stringsAsFactors = FALSE)$V1)

Lactobacillus_apis_orig_set <- read.table("../highqual_genomes_draft/Lactobacillus_apis.txt", stringsAsFactors = FALSE)$V1
Lactobacillus_apis_final_set <- c(Lactobacillus_apis_orig_set,
                                  "GCF_016102055.1", "GCA_003692845.1", "GCF_016101265.1", "GCF_016102085.1")
write.table(x = Lactobacillus_apis_final_set,
            file = "../highqual_genomes/Lactobacillus_apis.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

Lactobacillus_helsingborgensis_orig_set <- read.table("../highqual_genomes_draft/Lactobacillus_helsingborgensis.txt", stringsAsFactors = FALSE)$V1
Lactobacillus_helsingborgensis_final_set <- c(Lactobacillus_helsingborgensis_orig_set,
                                              "GCA_000761135.1", "GCA_003692925.1", "GCF_016100925.1")
write.table(x = Lactobacillus_helsingborgensis_final_set,
            file = "../highqual_genomes/Lactobacillus_helsingborgensis.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)


Lactobacillus_melliventris_orig_set <- read.table("../highqual_genomes_draft/Lactobacillus_melliventris.txt", stringsAsFactors = FALSE)$V1
Lactobacillus_melliventris_final_set <- c(Lactobacillus_melliventris_orig_set,
                                          "GCA_003692935.1", "GCF_016102065.1", "GCF_016102025.1",
                                          "GCF_016102045.1", "GCA_003693045.1")
write.table(x = Lactobacillus_melliventris_final_set,
            file = "../highqual_genomes/Lactobacillus_melliventris.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)


# These three genomes are closest to Lactobacillus kullabergensis, but are also very close to Lactobacillus kimbladii (only about 1 ANI closer).
# Decided to just exclude these genomes.
lactobacilli_dist_df["Lactobacillus sp|GCA_000760615.1", ]
lactobacilli_dist_df["Lactobacillus sp|GCF_016101625.1", ]
lactobacilli_dist_df["Lactobacillus sp|GCA_003693025.1", ]

# Similarly, these genomes are closer to L. kimbladii, but are ~92.4 ANI, so it's not clear. Decided to drop them as well.
lactobacilli_dist_df["Lactobacillus sp|GCF_016100885.1", ]
lactobacilli_dist_df["Lactobacillus sp|GCF_016100935.1", ]
lactobacilli_dist_df["Lactobacillus sp|GCF_016100975.1", ]


Lactobacillus_sp_to_exclude <- c("GCA_000760615.1", "GCF_016101625.1", "GCA_003693025.1",
                                 "GCF_016100885.1", "GCF_016100935.1", "GCF_016100975.1")


Lactobacillus_kullabergensis_final_set <- read.table("../highqual_genomes_draft/Lactobacillus_kullabergensis.txt", stringsAsFactors = FALSE)$V1
write.table(x = Lactobacillus_kullabergensis_final_set,
            file = "../highqual_genomes/Lactobacillus_kullabergensis.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

Lactobacillus_kimbladii_final_set <- read.table("../highqual_genomes_draft/Lactobacillus_kimbladii.txt", stringsAsFactors = FALSE)$V1
write.table(x = Lactobacillus_kimbladii_final_set,
            file = "../highqual_genomes/Lactobacillus_kimbladii.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)



#### Check that all Lactobacillus sp. accounted for.
Lactobacillus_sp_orig_set <- read.table("../highqual_genomes_draft/Lactobacillus_sp.txt", stringsAsFactors = FALSE)$V

Lactobacillus_sp_unaccounted <- Lactobacillus_sp_orig_set[-which(Lactobacillus_sp_orig_set %in% Lactobacillus_sp_to_exclude)]

Lactobacillus_sp_recoded <- c("GCA_003692935.1", "GCF_016102065.1", "GCF_016102025.1", "GCF_016102045.1",
                              "GCA_003693045.1", "GCA_000761135.1", "GCA_003692925.1", "GCF_016100925.1",
                              "GCF_016102055.1", "GCA_003692845.1", "GCF_016101265.1", "GCF_016102085.1",
                              "GCF_016100965.1", "GCF_016100985.1", "GCF_016101055.1", "GCF_016102125.1",
                              "GCF_016101025.1", "GCF_016101075.1", "GCF_016101045.1", "GCF_016101275.1")
Lactobacillus_sp_unaccounted <- Lactobacillus_sp_unaccounted[-which(Lactobacillus_sp_unaccounted %in% Lactobacillus_sp_recoded)]
Lactobacillus_sp_unaccounted



### Serratia_marcescens

Serratia_marcescens <- fastani_by_species$Serratia_marcescens
Serratia_marcescens <- Serratia_marcescens[-which(duplicated(Serratia_marcescens$G_combined)), ]
Serratia_marcescens <- Serratia_marcescens[-which(Serratia_marcescens$id1 == Serratia_marcescens$id2), ]
ggplot(data = Serratia_marcescens, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

Serratia_marcescens_dist_df <- generate_fastani_dist_df(Serratia_marcescens)
Serratia_marcescens_fastani_phylo <- ape::njs(as.dist(100 - Serratia_marcescens_dist_df))
plot(Serratia_marcescens_fastani_phylo)

Serratia_marcescens_final_set <- read.table("../highqual_genomes_draft/Serratia_marcescens.txt", stringsAsFactors = FALSE)$V1
write.table(x = Serratia_marcescens_final_set,
            file = "../highqual_genomes/Serratia_marcescens.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)


### Snodgrassella_alvi

Snodgrassella_alvi <- fastani_by_species$Snodgrassella_alvi
Snodgrassella_alvi <- Snodgrassella_alvi[-which(duplicated(Snodgrassella_alvi$G_combined)), ]
Snodgrassella_alvi <- Snodgrassella_alvi[-which(Snodgrassella_alvi$id1 == Snodgrassella_alvi$id2), ]
ggplot(data = Snodgrassella_alvi, aes(y = sp_combined, x = ani)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(cex=1, groupOnX = FALSE) +
  xlab("ANI") +
  ylab("")

Snodgrassella_alvi_dist_df <- generate_fastani_dist_df(Snodgrassella_alvi)
Snodgrassella_alvi_fastani_phylo <- ape::njs(as.dist(100 - Snodgrassella_alvi_dist_df))
plot(Snodgrassella_alvi_fastani_phylo)

Snodgrassella_alvi_final_set <- read.table("../highqual_genomes_draft/Snodgrassella_alvi.txt", stringsAsFactors = FALSE)$V1
write.table(x = Snodgrassella_alvi_final_set,
            file = "../highqual_genomes/Snodgrassella_alvi.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)


#### Also looked at dashing output, but it turned out to be a bit noisey and made it hard to cluster the genomes actually!
#### This may have been because some of the genomes had super high distances, approaching 1!

dashing_out <- read.table("dashing/distance_matrix.txt", header = FALSE, sep = "\t", row.names = 1)

rownames(dashing_out) <- gsub("../highqual_genomes/", "", rownames(dashing_out))
rownames(dashing_out) <- gsub(".fna", "", rownames(dashing_out))
rownames(dashing_out) <- gsub("/", "|", rownames(dashing_out))

colnames(dashing_out) <- rownames(dashing_out)

dashing_out[dashing_out == "-"] <- NA

dashing_out <- apply(dashing_out, 2, as.numeric)

diag(dashing_out) <- 0

dashing_out <- as.matrix(dashing_out)

dashing_out[lower.tri(dashing_out)] <- t(dashing_out)[lower.tri(dashing_out)]

dashing_out_similarity <- 1 - dashing_out

unique_23_col <- c("#4f3283",
                  "#92b540",
                  "#6c71d8",
                  "#cd9c2e",
                  "#c069c9",
                  "#63c671",
                  "#a74090",
                  "#45bc8d",
                  "#d94c53",
                  "#36dee6",
                  "#ca6333",
                  "#5694de",
                  "#c6944b",
                  "#7a70bd",
                  "#528739",
                  "#d48dd5",
                  "#93913d",
                  "#89285c",
                  "#d87c5c",
                  "#da6295",
                  "#893518",
                  "#b9475f",
                  "#b2464b")

species_names <- gsub("\\|.*$", "", colnames(dashing_out_similarity))

names(unique_23_col) <- sort(unique(species_names))

species_names_annot <- HeatmapAnnotation(species = species_names,
                                   col = list(species = unique_23_col))

Heatmap(dashing_out_similarity, show_row_names = FALSE, show_column_names = FALSE, top_annotation = species_names_annot)


dashing_out_dist <- as.dist(dashing_out)

dashing_out_phylo <- ape::nj(dashing_out_dist)

plot(dashing_out_phylo, "unrooted", main="NJ")

dashing_out_hclust <- hclust(d = dist(dashing_out), method = "complete")

