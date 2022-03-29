### Identify set of highquality genome accessions to use per species, based on CheckM output.

rm(list = ls(all.names = TRUE))

setwd("/data1/gdouglas/projects/honey_bee/ref_genomes/checkm_output/")

orig_checkm_output <- read.table("Carrie_orig_genome_checkm_output.tsv",
                                 header = TRUE, sep ="\t", stringsAsFactors = FALSE)

new_checkm_output <- read.table("2021_12_19_checkm_additional_downloads.tsv",
                                header = TRUE, sep ="\t", stringsAsFactors = FALSE)
new_checkm_output <- new_checkm_output[, -which(colnames(new_checkm_output) == "Marker.lineage.id")]

new_checkm_output_Lactobacillus_sp <- read.table("2022_01_05_checkm_additional_Lactobacillus_sp.txt",
                                header = TRUE, sep ="\t", stringsAsFactors = FALSE)
new_checkm_output_Lactobacillus_sp <- new_checkm_output_Lactobacillus_sp[, -which(colnames(new_checkm_output_Lactobacillus_sp) == "Marker.Lineage.ID")]


colnames(orig_checkm_output) <- colnames(new_checkm_output)
colnames(new_checkm_output_Lactobacillus_sp) <- colnames(new_checkm_output)

checkm_output <- rbind(orig_checkm_output, new_checkm_output)
checkm_output <- rbind(checkm_output, new_checkm_output_Lactobacillus_sp)

checkm_output$accession <- gsub("^GCA_", "GCAXXXXX", checkm_output$Bin)
checkm_output$accession <- gsub("^GCF_", "GCFXXXXX", checkm_output$accession)
checkm_output$accession <- gsub("_.*$", "", checkm_output$accession)
checkm_output$accession <- gsub("XXXXX", "_", checkm_output$accession)


# Remove duplicated id ("GCF_900094785.1")
checkm_output <- checkm_output[-which(duplicated(checkm_output$accession)), ]


# Check that no accessions duplicated between RefSeq / Genbank (i.e., same id except GCF/GCA swapped)
checkm_output$accession_GCA_only <- gsub("GCF_", "GCA_", checkm_output$accession)
which(duplicated(checkm_output$accession_GCA_only))



rownames(checkm_output) <- checkm_output$accession_GCA_only

accessions_by_species <- list()
species_names <- gsub(".tsv$", "", list.files("../adding_new_microbiota_genomes/accessions_to_process/"))
accession_filepaths <- list.files("../adding_new_microbiota_genomes/accessions_to_process/", full.names = TRUE)
for (i in 1:length(accession_filepaths)) {
  accessions_by_species[[species_names[i]]] <- gsub("GCF_", "GCA_", read.table(accession_filepaths[i],
                                                          stringsAsFactors = FALSE, sep = "\t", header = TRUE)$accession)
}

checkm_output_highqual <- checkm_output[which(checkm_output$Completeness > 95 & checkm_output$Contamination < 2), ]
genomes_remaining <- data.frame(matrix(NA, nrow = length(species_names), ncol = 3))
rownames(genomes_remaining) <- species_names
colnames(genomes_remaining) <- c("total", "in_checkm", "highqual")


for (s in species_names) {
  species_accessions <- accessions_by_species[[s]]
  genomes_remaining[s, "total"] <- length(species_accessions)
  genomes_remaining[s, "in_checkm"] <- length(which(species_accessions %in% rownames(checkm_output)))
  genomes_remaining[s, "highqual"] <- length(which(species_accessions %in% rownames(checkm_output_highqual)))
  
  highqual_species_accessions <- species_accessions[which(species_accessions %in% rownames(checkm_output_highqual))]
  
  
  outfile <- paste0("/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes_draft/", s, ".txt")
  write.table(x = checkm_output_highqual[highqual_species_accessions, "accession"],
              file = outfile,
              quote = FALSE, col.names = FALSE, row.names = FALSE)
}

