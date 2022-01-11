rm(list = ls(all.names = TRUE))

library(xlsx)

# Parse accessions for all Apis mellifera microbiota genomes (for specified species).
# Take the union of genomes available through IMG / NCBI (but mark both ids - NA if missing).
# Check against the Ellegaard 2020 genome database and figure out if any other species / genomes should be added that didn't come up in my search.

# Read in Ellegaard 2020 Data S1.

sheet_info <- read.table("/data1/gdouglas/projects/honey_bee/ref_genomes/adding_new_microbiota_genomes/Ellegaard_2020_Data_S1/sheet_descript.txt", stringsAsFactors = FALSE, sep = "\t")$V1
Ellegaard_db <- list()
for (i in 1:length(sheet_info)) {
  Ellegaard_db[[sheet_info[[i]]]] <- data.frame(read.xlsx("/data1/gdouglas/projects/honey_bee/ref_genomes/adding_new_microbiota_genomes/Ellegaard_2020_Data_S1/1-s2.0-S0960982220305868-mmc2.xlsx",
                                                          sheetIndex = i), stringsAsFactors = FALE)
}


# Apibacter_adventoris

Apibacter_adventoris_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Apibacter_adventoris.txt",
                                       sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Apibacter_adventoris_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Apibacter_adventoris.csv",
                                        sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_db$`(F) Apibacter spp.`



# Apilactobacillus_apinorum

# Only one genome in NCBI for this species (same as one used in Ellegaard 2020)
# Added manually
Apilactobacillus_apinorum_out <- data.frame(accession = "GCA_001281175.1",
                                              NCBI_download = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/281/175/GCA_001281175.1_ASM128117v1/",
                                              databae = "NCBI")

write.table(x = Apilactobacillus_apinorum_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Apilactobacillus_apinorum.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Apilactobacillus_kunkeei

Apilactobacillus_kunkeei_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Apilactobacillus_kunkeei.txt",
                                       sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Apilactobacillus_kunkeei_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Apilactobacillus_kunkeei.csv",
                                        sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Apilactobacillus_kunkeei_NCBI_Amellifera <- Apilactobacillus_kunkeei_NCBI[which(Apilactobacillus_kunkeei_NCBI$Host == "Apis mellifera"), ]

Ellegaard_db$`(K) Lactobacillus spp.`

Apilactobacillus_kunkeei_out <- Apilactobacillus_kunkeei_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Apilactobacillus_kunkeei_out) <- c("accession", "NCBI_download")
Apilactobacillus_kunkeei_out$database <- "NCBI"

write.table(x = Apilactobacillus_kunkeei_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Apilactobacillus_kunkeei.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Bartonella_apis

Bartonella_apis_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Bartonella_apis.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Bartonella_apis_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Bartonella_apis.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_db$`(H) Bartonella apis`

Bartonella_apis_IMG_in_Ellegaard <- Bartonella_apis_IMG[which(Bartonella_apis_IMG$IMG.Genome.ID %in% Ellegaard_db$`(H) Bartonella apis`$IMG_acc), ]

Bartonella_apis_IMG_in_Ellegaard_in_NCBI <- Bartonella_apis_IMG_in_Ellegaard[which(Bartonella_apis_IMG_in_Ellegaard$NCBI.Assembly.Accession %in% Bartonella_apis_NCBI$Assembly), ] 


Bartonella_apis_NCBI_Amellifera <- Bartonella_apis_NCBI[which(Bartonella_apis_NCBI$Host == "Apis mellifera"), ]

Bartonella_apis_out <- Bartonella_apis_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Bartonella_apis_out) <- c("accession", "NCBI_download")
Bartonella_apis_out$database <- "NCBI"

write.table(x = Bartonella_apis_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Bartonella_apis.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Bifidobacterium_asteroides

Bifidobacterium_asteroides_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Bifidobacterium_asteroides.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Bifidobacterium_asteroides_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Bifidobacterium_asteroides.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_db$`(A) Bifidobacterium spp.`

Bifidobacterium_asteroides_IMG_in_Ellegaard <- Bifidobacterium_asteroides_IMG[which(Bifidobacterium_asteroides_IMG$IMG.Genome.ID %in% Ellegaard_db$`(A) Bifidobacterium spp.`$IMG_acc), ]

Bifidobacterium_asteroides_IMG_in_Ellegaard_in_NCBI <- Bifidobacterium_asteroides_IMG_in_Ellegaard[which(Bifidobacterium_asteroides_IMG_in_Ellegaard$NCBI.Assembly.Accession %in% Bifidobacterium_asteroides_NCBI$Assembly), ] 

Bifidobacterium_asteroides_IMG_in_Ellegaard_NOT_in_NCBI <- Bifidobacterium_asteroides_IMG_in_Ellegaard[which(! Bifidobacterium_asteroides_IMG_in_Ellegaard$NCBI.Assembly.Accession %in% Bifidobacterium_asteroides_NCBI$Assembly), ] 

Ellegaard_db$`(A) Bifidobacterium spp.`[which(Ellegaard_db$`(A) Bifidobacterium spp.`$IMG_acc %in% Bifidobacterium_asteroides_IMG_in_Ellegaard_NOT_in_NCBI$taxon_oid), ]

Bifidobacterium_asteroides_NCBI_Amellifera <- Bifidobacterium_asteroides_NCBI[which(Bifidobacterium_asteroides_NCBI$Host == "Apis mellifera"), ]

# Some were present in Carrie's output that she found evidence that they were from Apis mellifera
Bifidobacterium_asteroides_NCBI_additional <- Bifidobacterium_asteroides_NCBI[which(Bifidobacterium_asteroides_NCBI$Assembly %in% c("GCA_002715865.1", "GCA_002846895.1")), ]

Bifidobacterium_asteroides_NCBI_Amellifera <- rbind(Bifidobacterium_asteroides_NCBI_Amellifera, Bifidobacterium_asteroides_NCBI_additional)

Bifidobacterium_asteroides_out <- Bifidobacterium_asteroides_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Bifidobacterium_asteroides_out) <- c("accession", "NCBI_download")
Bifidobacterium_asteroides_out$database <- "NCBI"

# Note that I add this genome manually which was listed as "Bifidobacterium asteroides PRL2011" on NCBI and which was present in the Ellegaard database.
Bifidobacterium_asteroides_out[15, ] <- c("GCA_000304215.1", "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/304/215/GCA_000304215.1_ASM30421v1", "NCBI")

# Also added these two genomes whichI think are likely mislabelled, but I need to check:
Bifidobacterium_asteroides_out[16, ] <- c("GCA_007559155.1", "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/559/155/GCA_007559155.1_ASM755915v1", "NCBI")
Bifidobacterium_asteroides_out[17, ] <- c("GCA_007559275.1", "ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/559/275/GCA_007559275.1_ASM755927v1", "NCBI")

write.table(x = Bifidobacterium_asteroides_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Bifidobacterium_asteroides.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


# Bifidobacterium_coryneforme

Bifidobacterium_coryneforme_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Bifidobacterium_coryneforme.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Bifidobacterium_coryneforme_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Bifidobacterium_coryneforme.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")



Ellegaard_db_coryneforme <- Ellegaard_db$`(A) Bifidobacterium spp.`[grep("coryneforme", Ellegaard_db$`(A) Bifidobacterium spp.`$Species), ]

Bifidobacterium_coryneforme_IMG_in_Ellegaard <- Bifidobacterium_coryneforme_IMG[which(Bifidobacterium_coryneforme_IMG$IMG.Genome.ID %in% Ellegaard_db$`(H) Bartonella apis`$IMG_acc), ]

Bifidobacterium_coryneforme_IMG_in_Ellegaard_in_NCBI <- Bifidobacterium_coryneforme_IMG_in_Ellegaard[which(Bifidobacterium_coryneforme_IMG_in_Ellegaard$NCBI.Assembly.Accession %in% Bifidobacterium_coryneforme_NCBI$Assembly), ] 


# After some searching was able to determine that strain LMG18911 was from Apis mellifera based on linking ids to Ellegaard 2020, but DSM 20216 is ambiguous as the paper where it was isolated in 1969 is inaccessible and they looked at multiple honey bee species.
Bifidobacterium_coryneforme_NCBI_Amellifera <- Bifidobacterium_coryneforme_NCBI[c(1, 3), ]

Bifidobacterium_coryneforme_out <- Bifidobacterium_coryneforme_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Bifidobacterium_coryneforme_out) <- c("accession", "NCBI_download")
Bifidobacterium_coryneforme_out$database <- "NCBI"

write.table(x = Bifidobacterium_coryneforme_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Bifidobacterium_coryneforme.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Bifidobacterium_indicum

Bifidobacterium_indicum_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Bifidobacterium_indicum.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Bifidobacterium_indicum_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Bifidobacterium_indicum.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")


Bifidobacterium_indicum_NCBI_Amellifera <- Bifidobacterium_indicum_NCBI[which(Bifidobacterium_indicum_NCBI$Host == "Apis mellifera"), ]

Bifidobacterium_indicum_out <- Bifidobacterium_indicum_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Bifidobacterium_indicum_out) <- c("accession", "NCBI_download")
Bifidobacterium_indicum_out$database <- "NCBI"

write.table(x = Bifidobacterium_indicum_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Bifidobacterium_indicum.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


# Bombella_apis

Bombella_apis_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Bombella_apis.csv",
                                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Bombella_apis_NCBI_Amellifera <- Bombella_apis_NCBI[which(Bombella_apis_NCBI$Host == "Apis mellifera"), ]

Bombella_apis_out <- Bombella_apis_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Bombella_apis_out) <- c("accession", "NCBI_download")
Bombella_apis_out$database <- "NCBI"

write.table(x = Bombella_apis_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Bombella_apis.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


# Bombella_sp

Bombella_sp_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Bombella_sp.csv",
                                 sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Bombella_sp_NCBI_Amellifera <- Bombella_sp_NCBI[which(Bombella_sp_NCBI$Host == "Apis mellifera"), ]

Bombella_sp_out <- Bombella_sp_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Bombella_sp_out) <- c("accession", "NCBI_download")
Bombella_sp_out$database <- "NCBI"

# needed to add some genbank links by hand:
Bombella_sp_out[1, "NCBI_download"] <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/734/205/GCA_009734205.1_ASM973420v1"
Bombella_sp_out[3, "NCBI_download"] <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/725/745/GCA_009725745.1_ASM972574v1"
Bombella_sp_out[4, "NCBI_download"] <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/725/805/GCA_009725805.1_ASM972580v1"

write.table(x = Bombella_sp_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Bombella_sp.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


# Bombilactobacillus_mellifer

Bombilactobacillus_mellifer_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Bombilactobacillus_mellifer.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

# Only one accession in this case.
Bombilactobacillus_mellifer_NCBI <- "GCA_000970795.1"


Ellegaard_db$`(C) Lactobacillus spp., phylotype Firm4`

# Write out this single accession, which is the only one by all three approaches.
Bombilactobacillus_mellifer_out <- data.frame(accession = "GCA_000970795.1",
                                              NCBI_download = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/970/795/GCA_000970795.1_ASM97079v1/",
                                              databae = "NCBI")

write.table(x = Bombilactobacillus_mellifer_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Bombilactobacillus_mellifer.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Bombilactobacillus_mellis

Bombilactobacillus_mellis_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Bombilactobacillus_mellis.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Bombilactobacillus_mellis_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Bombilactobacillus_mellis.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_db$`(C) Lactobacillus spp., phylotype Firm4`

Bombilactobacillus_mellis_IMG_in_Ellegaard <- Bombilactobacillus_mellis_IMG[which(Bombilactobacillus_mellis_IMG$IMG.Genome.ID %in% Ellegaard_db$`(C) Lactobacillus spp., phylotype Firm4`$IMG_acc), ]

Bombilactobacillus_mellis_IMG_in_Ellegaard_in_NCBI <- Bombilactobacillus_mellis_IMG_in_Ellegaard[which(Bombilactobacillus_mellis_IMG_in_Ellegaard$NCBI.Assembly.Accession %in% Bombilactobacillus_mellis_NCBI$Assembly), ] 

Ellegaard_db$`(C) Lactobacillus spp., phylotype Firm4`[which(Ellegaard_db$`(C) Lactobacillus spp., phylotype Firm4`$IMG_acc %in% Bombilactobacillus_mellis_IMG$taxon_oid), ]

Bombilactobacillus_mellis_NCBI_Amellifera <- Bombilactobacillus_mellis_NCBI[which(Bombilactobacillus_mellis_NCBI$Host == "Apis mellifera"), ]

Bombilactobacillus_mellis_out <- Bombilactobacillus_mellis_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Bombilactobacillus_mellis_out) <- c("accession", "NCBI_download")
Bombilactobacillus_mellis_out$database <- "NCBI"

write.table(x = Bombilactobacillus_mellis_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Bombilactobacillus_mellis.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


#Commensalibacter_sp

Commensalibacter_sp_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Commensalibacter_sp.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Commensalibacter_sp_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Commensalibacter_sp.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_db$`(J) Commensalibacter spp.`

Commensalibacter_sp_IMG_in_Ellegaard <- Commensalibacter_sp_IMG[which(Commensalibacter_sp_IMG$IMG.Genome.ID %in% Ellegaard_db$`(H) Bartonella apis`$IMG_acc), ]

Commensalibacter_sp_IMG_in_Ellegaard_in_NCBI <- Commensalibacter_sp_IMG_in_Ellegaard[which(Commensalibacter_sp_IMG_in_Ellegaard$NCBI.Assembly.Accession %in% Commensalibacter_sp_NCBI$Assembly), ] 



Commensalibacter_sp_NCBI_Amellifera <- Commensalibacter_sp_NCBI[which(Commensalibacter_sp_NCBI$Host == "Apis mellifera"), ]

Commensalibacter_sp_out <- Commensalibacter_sp_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Commensalibacter_sp_out) <- c("accession", "NCBI_download")
Commensalibacter_sp_out$database <- "NCBI"

write.table(x = Commensalibacter_sp_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Commensalibacter_sp.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Frischella_perrara

Frischella_perrara_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Frischella_perrara.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Frischella_perrara_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Frischella_perrara.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_db$`(G) Frischella perrrara`

Frischella_perrara_IMG_in_Ellegaard <- Frischella_perrara_IMG[which(Frischella_perrara_IMG$IMG.Genome.ID %in% Ellegaard_db$`(H) Bartonella apis`$IMG_acc), ]

Frischella_perrara_IMG_in_Ellegaard_in_NCBI <- Frischella_perrara_IMG_in_Ellegaard[which(Frischella_perrara_IMG_in_Ellegaard$NCBI.Assembly.Accession %in% Frischella_perrara_NCBI$Assembly), ] 



Frischella_perrara_NCBI_Amellifera <- Frischella_perrara_NCBI[which(Frischella_perrara_NCBI$Host == "Apis mellifera"), ]

Frischella_perrara_out <- Frischella_perrara_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Frischella_perrara_out) <- c("accession", "NCBI_download")
Frischella_perrara_out$database <- "NCBI"

write.table(x = Frischella_perrara_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Frischella_perrara.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Fructobacillus_fructosus

Fructobacillus_fructosus_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Fructobacillus_fructosus.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Fructobacillus_fructosus_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Fructobacillus_fructosus.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_db$Fructobacillus_fructosus

# Non-applicable


# Gilliamella_apicola

Gilliamella_apicola_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Gilliamella_apicola.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Gilliamella_apicola_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Gilliamella_apicola.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_Gilliamella_apicola <- Ellegaard_db$`(D) Gilliamella spp.`[grep("Gilliamella apicola", Ellegaard_db$`(D) Gilliamella spp.`$Species), ]

Gilliamella_apicola_IMG_in_Ellegaard <- Gilliamella_apicola_IMG[which(Gilliamella_apicola_IMG$IMG.Genome.ID %in% Ellegaard_db$`(H) Bartonella apis`$IMG_acc), ]

Gilliamella_apicola_IMG_in_Ellegaard_in_NCBI <- Gilliamella_apicola_IMG_in_Ellegaard[which(Gilliamella_apicola_IMG_in_Ellegaard$NCBI.Assembly.Accession %in% Gilliamella_apicola_NCBI$Assembly), ] 



Gilliamella_apicola_NCBI_Amellifera <- Gilliamella_apicola_NCBI[which(Gilliamella_apicola_NCBI$Host == "Apis mellifera"), ]

Gilliamella_apicola_NCBI_additional <- Gilliamella_apicola_NCBI[which(Gilliamella_apicola_NCBI$Assembly == "GCA_000599985.1"),]

Gilliamella_apicola_NCBI_Amellifera <- rbind(Gilliamella_apicola_NCBI_Amellifera, Gilliamella_apicola_NCBI_additional)

Gilliamella_apicola_out <- Gilliamella_apicola_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Gilliamella_apicola_out) <- c("accession", "NCBI_download")
Gilliamella_apicola_out$database <- "NCBI"

Gilliamella_apicola_no_refseq_accessions <- Gilliamella_apicola_out[which(Gilliamella_apicola_out$NCBI_download == ""), "accession"]

Gilliamella_apicola_out[which(Gilliamella_apicola_out$accession == "GCA_001723875.1"), "NCBI_download"] <- Gilliamella_apicola_NCBI_Amellifera[which(Gilliamella_apicola_NCBI_Amellifera$Assembly == "GCA_001723875.1"), "GenBank.FTP"]
Gilliamella_apicola_out[which(Gilliamella_apicola_out$accession == "GCA_001726005.1"), "NCBI_download"] <- Gilliamella_apicola_NCBI_Amellifera[which(Gilliamella_apicola_NCBI_Amellifera$Assembly == "GCA_001726005.1"), "GenBank.FTP"]

write.table(x = Gilliamella_apicola_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Gilliamella_apicola.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)




# Gilliamella_apis

Gilliamella_apis_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Gilliamella_apis.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Gilliamella_apis_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Gilliamella_apis.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_Gilliamella_apis<- Ellegaard_db$`(D) Gilliamella spp.`[grep("Gilliamella apis", Ellegaard_db$`(D) Gilliamella spp.`$Species), ]

Ellegaard_Gilliamella_apis$Strain_name <- gsub(" $", "", Ellegaard_Gilliamella_apis$Strain_name)


Gilliamella_apis_NCBI_Amellifera <- Gilliamella_apis_NCBI[which(Gilliamella_apis_NCBI$Host == "Apis mellifera"), ]

which(! Ellegaard_Gilliamella_apis$Strain_name %in% Gilliamella_apis_NCBI_Amellifera$Strain)

Gilliamella_apis_out <- Gilliamella_apis_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Gilliamella_apis_out) <- c("accession", "NCBI_download")
Gilliamella_apis_out$database <- "NCBI"

write.table(x = Gilliamella_apis_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Gilliamella_apis.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Hafnia_alvei

Hafnia_alvei_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Hafnia_alvei.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Hafnia_alvei_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Hafnia_alvei.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

# No Apis mellifera found



# Lactobacillus_apis

Lactobacillus_apis_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Lactobacillus_apis.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Lactobacillus_apis_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Lactobacillus_apis.csv",
                                      sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_db$`(B) Lactobacillus spp., phylotype Firm5`
# Strains in the database
# ESL0185
# R-53131/LMG 26964
# ESL0353/Hma11
# ESL0263

Lactobacillus_apis_IMG_in_Ellegaard <- Lactobacillus_apis_IMG[which(Lactobacillus_apis_IMG$IMG.Genome.ID %in% Ellegaard_db$`(H) Bartonella apis`$IMG_acc), ]

Lactobacillus_apis_IMG_in_Ellegaard_in_NCBI <- Lactobacillus_apis_IMG_in_Ellegaard[which(Lactobacillus_apis_IMG_in_Ellegaard$NCBI.Assembly.Accession %in% Lactobacillus_apis_NCBI$Assembly), ] 


# Based on Carrie's notes I know that GCA_900094785.1 is associated with Apis mellifera too. So can keep all except for one that I'm unsure about.

Lactobacillus_apis_NCBI_Amellifera <- Lactobacillus_apis_NCBI[-2, ]

Lactobacillus_apis_out <- Lactobacillus_apis_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Lactobacillus_apis_out) <- c("accession", "NCBI_download")
Lactobacillus_apis_out$database <- "NCBI"

write.table(x = Lactobacillus_apis_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Lactobacillus_apis.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Lactobacillus_helsingborgensis

Lactobacillus_helsingborgensis_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Lactobacillus_helsingborgensis.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Lactobacillus_helsingborgensis_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Lactobacillus_helsingborgensis.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_db$`(B) Lactobacillus spp., phylotype Firm5`

# Strains:
# ESL0183
# ESL0354/Bma5
# wkB8

Lactobacillus_helsingborgensis_NCBI_Amellifera <- Lactobacillus_helsingborgensis_NCBI[which(Lactobacillus_helsingborgensis_NCBI$Host == "Apis mellifera"), ]

Lactobacillus_helsingborgensis_out <- Lactobacillus_helsingborgensis_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Lactobacillus_helsingborgensis_out) <- c("accession", "NCBI_download")
Lactobacillus_helsingborgensis_out$database <- "NCBI"

write.table(x = Lactobacillus_helsingborgensis_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Lactobacillus_helsingborgensis.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)




# Lactobacillus_kimbladii

Lactobacillus_kimbladii_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Lactobacillus_kimbladii.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Lactobacillus_kimbladii_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Lactobacillus_kimbladii.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_db$`(B) Lactobacillus spp., phylotype Firm5`
# Only strain:
# 2627854222    JXLH00000000.1  Ga0072403           JF75     ESL0352             Hma2     Firm5 Firm5-4        Lactobacillus kimbladii    Apis mellifera 


Lactobacillus_kimbladii_NCBI_Amellifera <- Lactobacillus_kimbladii_NCBI[which(Lactobacillus_kimbladii_NCBI$Host == "Apis mellifera"), ]

Lactobacillus_kimbladii_out <- Lactobacillus_kimbladii_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Lactobacillus_kimbladii_out) <- c("accession", "NCBI_download")
Lactobacillus_kimbladii_out$database <- "NCBI"

write.table(x = Lactobacillus_kimbladii_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Lactobacillus_kimbladii.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Lactobacillus_kullabergensis

Lactobacillus_kullabergensis_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Lactobacillus_kullabergensis.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Lactobacillus_kullabergensis_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Lactobacillus_kullabergensis.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_db$`(B) Lactobacillus spp., phylotype Firm5`
# Four strains:
# 14 2636415759    JXBY00000000.1  Ga0070887           JF76     ESL0351            Biut2     Firm5 Firm5-4   Lactobacillus kullabergensis    Apis mellifera                       [S2]
# 15         NA      JRJB00000000   LACWKB10           <NA>       wkB10             <NA>     Firm5 Firm5-4   Lactobacillus kullabergensis    Apis mellifera                       [S7]
# 16 2758568514      REHN00000000  Ga0226839           <NA>     ESL0261             <NA>     Firm5 Firm5-4   Lactobacillus kullabergensis    Apis mellifera                       [S6]
# 17 2684622911          CP029477  Ga0133564           <NA>     ESL0186             <NA>     Firm5 Firm5-4  Lactobacillus kulllabergensis    Apis mellifera                       [S1]



Lactobacillus_kullabergensis_NCBI_Amellifera <- Lactobacillus_kullabergensis_NCBI[which(Lactobacillus_kullabergensis_NCBI$Host == "Apis mellifera"), ]

Lactobacillus_kullabergensis_out <- Lactobacillus_kullabergensis_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Lactobacillus_kullabergensis_out) <- c("accession", "NCBI_download")
Lactobacillus_kullabergensis_out$database <- "NCBI"

write.table(x = Lactobacillus_kullabergensis_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Lactobacillus_kullabergensis.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Lactobacillus_melliventris

Lactobacillus_melliventris_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Lactobacillus_melliventris.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Lactobacillus_melliventris_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Lactobacillus_melliventris.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_db$`(B) Lactobacillus spp., phylotype Firm5`
# 8  2684622913      QGLG00000000  Ga0133562           <NA>     ESL0184             <NA>     Firm5 Firm5-3     Lactobacillus melliventris    Apis mellifera                       [S1]
# 9  2639762991    JXLI00000000.1  Ga0072404           JF74     ESL0350             Hma8     Firm5 Firm5-3     Lactobacillus melliventris    Apis mellifera                       [S2]
# 10 2758568515      REHP00000000  Ga0225908           <NA>     ESL0259             <NA>     Firm5 Firm5-3     Lactobacillus melliventris    Apis mellifera                       [S6]
# 11 2758568513      REHO00000000  Ga0226840           <NA>     ESL0260             <NA>     Firm5 Firm5-3     Lactobacillus melliventris    Apis mellifera                       [S6]
# 12 2758568558               NA   Ga0227359           <NA>     ESL0393             <NA>     Firm5 Firm5-3     Lactobacillus melliventris    Apis mellifera


Lactobacillus_melliventris_NCBI_Amellifera <- Lactobacillus_melliventris_NCBI[which(Lactobacillus_melliventris_NCBI$Host == "Apis mellifera"), ]

Lactobacillus_melliventris_out <- Lactobacillus_melliventris_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Lactobacillus_melliventris_out) <- c("accession", "NCBI_download")
Lactobacillus_melliventris_out$database <- "NCBI"

write.table(x = Lactobacillus_melliventris_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Lactobacillus_melliventris.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


# Lactobacillus sp.

Lactobacillus_sp_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/ref_genomes/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Lactobacillus_sp.csv",
                                              sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Lactobacillus_sp_NCBI_Amellifera <- Lactobacillus_sp_NCBI[which(Lactobacillus_sp_NCBI$Host == "Apis mellifera"), ]

Lactobacillus_sp_out <- Lactobacillus_sp_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Lactobacillus_sp_out) <- c("accession", "NCBI_download")
Lactobacillus_sp_out$database <- "NCBI"

write.table(x = Lactobacillus_sp_out,
            file = "/data1/gdouglas/projects/honey_bee/ref_genomes/adding_new_microbiota_genomes/accessions_to_process/Lactobacillus_sp.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Parasaccharibacter_apium

Parasaccharibacter_apium_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Parasaccharibacter_apium.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Parasaccharibacter_apium_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Parasaccharibacter_apium.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

# Even the one not listed as Apis mellifera was isolated from Apis mellifera, which I determined by looking at the associated bioproject.
Parasaccharibacter_apium_NCBI_Amellifera <- Parasaccharibacter_apium_NCBI

Parasaccharibacter_apium_out <- Parasaccharibacter_apium_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Parasaccharibacter_apium_out) <- c("accession", "NCBI_download")
Parasaccharibacter_apium_out$database <- "NCBI"

write.table(x = Parasaccharibacter_apium_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Parasaccharibacter_apium.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Saccharibacter_sp

Saccharibacter_sp_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Saccharibacter_sp.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Saccharibacter_sp_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Saccharibacter_sp.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

# One genome with unlisted Host I found out through searching bioproject that it is from the same isolate from Apis mellifera as M18
# Veress, A., Wilk, T., Kiss, J., Olasz, F., and Papp, P.P. "Draft genome sequences of Saccharibacter sp. strains 3.A.1 and M18 isolated from honey and a honey bee (Apis mellifera) stomach." Genome Announc. (2017) 5:e00744-17.
Saccharibacter_sp_NCBI_Amellifera <- Saccharibacter_sp_NCBI[which(Saccharibacter_sp_NCBI$Strain %in% c("M18", "3.A.1")), ]

Saccharibacter_sp_out <- Saccharibacter_sp_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Saccharibacter_sp_out) <- c("accession", "NCBI_download")
Saccharibacter_sp_out$database <- "NCBI"

write.table(x = Saccharibacter_sp_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Saccharibacter_sp.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Serratia_marcescens

Serratia_marcescens_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Serratia_marcescens.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Serratia_marcescens_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Serratia_marcescens.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")


Serratia_marcescens_NCBI_Amellifera <- Serratia_marcescens_NCBI[which(Serratia_marcescens_NCBI$Host == "Apis mellifera"), ]

Serratia_marcescens_out <- Serratia_marcescens_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Serratia_marcescens_out) <- c("accession", "NCBI_download")
Serratia_marcescens_out$database <- "NCBI"

write.table(x = Serratia_marcescens_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Serratia_marcescens.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



# Snodgrassella_alvi

Snodgrassella_alvi_IMG <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/IMG_Snodgrassella_alvi.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Snodgrassella_alvi_NCBI <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_12_16_honey_bee_microbiota_info/NCBI_Snodgrassella_alvi.csv",
                       sep = ",", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

Ellegaard_Snodgrassella_alvi <- Ellegaard_db$`(E) Snodgrassella spp.`[which(Ellegaard_db$`(E) Snodgrassella spp.`$Species == "Snodgrassella alvi"), ]
Ellegaard_Snodgrassella_alvi$Strain_name <- gsub(" $", "", Ellegaard_Snodgrassella_alvi$Strain_name)
Ellegaard_Snodgrassella_alvi$Strain_name[which(! Ellegaard_Snodgrassella_alvi$Strain_name %in% Snodgrassella_alvi_NCBI$Strain)]
# Three strains not found in NCBI set:
# "ESL0323" "ESL0324" "ESL0304"
# 15 2758568599               NA Ga0227305     ESL0323 Snodgrassella Snod_1 Snodgrassella alvi Apis mellifera  This study                   Switzerland, Engel apiary
# 16 2758568598               NA Ga0227306     ESL0324 Snodgrassella Snod_1 Snodgrassella alvi Apis mellifera  This study                   Switzerland, Engel apiary
# 17 2758568600               NA Ga0227304     ESL0304 Snodgrassella Snod_1 Snodgrassella alvi  Apis mellifera This study                   Switzerland, Engel apiary


Snodgrassella_alvi_NCBI_Amellifera <- Snodgrassella_alvi_NCBI[which(Snodgrassella_alvi_NCBI$Host == "Apis mellifera"), ]

# Based on Carrie's annotations I realized that these two genomes should be included too:
Snodgrassella_alvi_NCBI_additional <- Snodgrassella_alvi_NCBI[which(Snodgrassella_alvi_NCBI$Assembly %in% c("GCA_000600005.1", "GCA_002406645.1")), ]

Snodgrassella_alvi_NCBI_Amellifera <- rbind(Snodgrassella_alvi_NCBI_Amellifera, Snodgrassella_alvi_NCBI_additional)

Snodgrassella_alvi_out <- Snodgrassella_alvi_NCBI_Amellifera[, c("Assembly", "RefSeq.FTP")]

colnames(Snodgrassella_alvi_out) <- c("accession", "NCBI_download")
Snodgrassella_alvi_out$database <- "NCBI"

write.table(x = Snodgrassella_alvi_out,
            file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/Snodgrassella_alvi.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)




