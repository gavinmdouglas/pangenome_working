### Process dbCAN2 output files (that correspond to the microbiota-wide panaroo pangenome, based on reference genomes only).
### There are eight output files as the full FASTA had to be split into eight parts to run on the dbCAN2 online portal.

rm(list = ls(all.names = TRUE))

library(seqinr)
library(stringr)

setwd("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_ref_only_annot/dbCAN2_raw_output_microbiota_panaroo_combined_protein_CDS")

input_CDS_seq_ids <- names(read.fasta("../../panaroo_out_ref_only/protein_CDS/microbiota_combined_protein_CDS.fasta.gz"))

dbCAN2_out <- list()
dbCAN2_out[["part1"]] <- read.table("dbCAN2_output_part1.txt.gz", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, comment.char = "", quote = "")
dbCAN2_out[["part2"]] <- read.table("dbCAN2_output_part2.txt.gz", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, comment.char = "", quote = "")
dbCAN2_out[["part3"]] <- read.table("dbCAN2_output_part3.txt.gz", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, comment.char = "", quote = "")
dbCAN2_out[["part4"]] <- read.table("dbCAN2_output_part4.txt.gz", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, comment.char = "", quote = "")
dbCAN2_out[["part5"]] <- read.table("dbCAN2_output_part5.txt.gz", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, comment.char = "", quote = "")
dbCAN2_out[["part6"]] <- read.table("dbCAN2_output_part6.txt.gz", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, comment.char = "", quote = "")
dbCAN2_out[["part7"]] <- read.table("dbCAN2_output_part7.txt.gz", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, comment.char = "", quote = "")
dbCAN2_out[["part8"]] <- read.table("dbCAN2_output_part8.txt.gz", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, comment.char = "", quote = "")

dbCAN2_all_hits <- do.call(what = rbind, args = dbCAN2_out)
rownames(dbCAN2_all_hits) <- gsub("^part\\d\\.", "", rownames(dbCAN2_all_hits))

# Percent of CDS sequences annotated with CAZy (based on at least one approach):
(nrow(dbCAN2_all_hits) / length(input_CDS_seq_ids)) * 100

# Only keep rows where value of # tools called >= 2
dbCAN2_all_hits_atleast2 <- dbCAN2_all_hits[which(dbCAN2_all_hits$X.ofTools >= 2), ]


# Clean up HMMER column
dbCAN2_all_hits_atleast2$HMMER <- gsub("\\(\\d*-\\d*\\)", "", dbCAN2_all_hits_atleast2$HMMER)
dbCAN2_all_hits_atleast2$HMMER <- gsub("GT2_Glycos_transf_2", "GT2", dbCAN2_all_hits_atleast2$HMMER)
dbCAN2_all_hits_atleast2$HMMER <- gsub("GT2_Glyco_tranf_2", "GT2", dbCAN2_all_hits_atleast2$HMMER)
dbCAN2_all_hits_atleast2$HMMER <- gsub("GT2_Glyco_trans_2", "GT2", dbCAN2_all_hits_atleast2$HMMER)

# Remove Singalp and X.ofTools columns 
dbCAN2_all_hits_atleast2 <- dbCAN2_all_hits_atleast2[, -which(colnames(dbCAN2_all_hits_atleast2) %in% c("Signalp", "X.ofTools"))]

# Get consensus (keep all calls by at least two tools)
dbCAN2_all_hits_atleast2$Consensus <- NA

for (i in 1:nrow(dbCAN2_all_hits_atleast2)) {
  
  HMMER_calls <- unique(str_split(dbCAN2_all_hits_atleast2[i, "HMMER"], "\\+")[[1]])
  Hotpep_calls <- unique(str_split(dbCAN2_all_hits_atleast2[i, "Hotpep"], "\\+")[[1]])
  DIAMOND_calls <- unique(str_split(dbCAN2_all_hits_atleast2[i, "DIAMOND"], "\\+")[[1]])
  
  individual_calls <- c(HMMER_calls, Hotpep_calls, DIAMOND_calls)
  
  if ("N" %in% individual_calls) { individual_calls <- individual_calls[-which(individual_calls == "N")] }
  
  individual_calls_tally <- table(individual_calls)
  
  singleton_cazy_calls <- names(individual_calls_tally[which(individual_calls_tally == 1)])
  multi_calls <- names(individual_calls_tally[which(individual_calls_tally >= 2)])
  
  # If any singletons are subfamilies then try removing the subfamily (_x) prefix just in case that leads to more consensuses at the gene family level.
  if ((length(grep("_", singleton_cazy_calls)) > 0) & (length(singleton_cazy_calls) > 1)) {
    singleton_cazy_calls <- gsub("_.*$", "", singleton_cazy_calls)
    singleton_cazy_calls_retallied <- table(singleton_cazy_calls)
    multi_calls <- c(multi_calls, names(singleton_cazy_calls_retallied[which(singleton_cazy_calls_retallied >= 2)]))
  }
  
  dbCAN2_all_hits_atleast2[i, "Consensus"] <- paste(sort(multi_calls), collapse = ",")

}


# Remove remaining entries where no consensus could be found (despite the above fixing and the X.ofTools column saying at least 2)
dbCAN2_all_hits_atleast2 <- dbCAN2_all_hits_atleast2[-which(dbCAN2_all_hits_atleast2$Consensus == ""), ]

# Get percent that could be called by at least 2 tools (consistently):
(nrow(dbCAN2_all_hits_atleast2) / length(input_CDS_seq_ids)) * 100

# Set all values of "N" to be NA
dbCAN2_all_hits_atleast2$HMMER[which(dbCAN2_all_hits_atleast2$HMMER == "N")] <- NA
dbCAN2_all_hits_atleast2$Hotpep[which(dbCAN2_all_hits_atleast2$Hotpep == "N")] <- NA
dbCAN2_all_hits_atleast2$DIAMOND[which(dbCAN2_all_hits_atleast2$DIAMOND == "N")] <- NA

dbCAN2_all_hits_atleast2 <- dbCAN2_all_hits_atleast2[sort(rownames(dbCAN2_all_hits_atleast2)), ]

write.table(x = dbCAN2_all_hits_atleast2,
            file = "../dbCAN2_microbiota_panaroo_combined_protein_CDS_consensus.tsv",
            row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
