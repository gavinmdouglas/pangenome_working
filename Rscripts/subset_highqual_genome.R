### Restrict ids in genome accession files to those marked as high quality based on checkm.

rm(list = ls(all.names = TRUE))

setwd("/Users/Gavin/Google_Drive/postdoc/honey_bee_pangenome/data/microbiota_genomes/")

library(data.table)

checkm_output <- fread("checkm_output.fixed.tsv")

checkm_output_highqual <- checkm_output[completeness >= 97 & contamination <= 2]
checkm_output_lowqual <- checkm_output[completeness < 97 | contamination > 2]

checkm_output_highqual$accession <- gsub("_genomic", "", checkm_output_highqual$bin)

orig_Bifidobacterium_ids <- read.table("all_genome_ids/Bifidobacterium_ids.txt", stringsAsFactors = FALSE)$V1
highqual_Bifidobacterium_ids <- orig_Bifidobacterium_ids[which(orig_Bifidobacterium_ids %in% checkm_output_highqual$accession)]
lowqual_Bifidobacterium_ids <- orig_Bifidobacterium_ids[-which(orig_Bifidobacterium_ids %in% checkm_output_highqual$accession)]
write.table(x = highqual_Bifidobacterium_ids, file = "highqual_genome_ids/Bifidobacterium_ids.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

orig_Firm4_ids <- read.table("all_genome_ids/Firm4_ids.txt", stringsAsFactors = FALSE)$V1
highqual_Firm4_ids <- orig_Firm4_ids[which(orig_Firm4_ids %in% checkm_output_highqual$accession)]
lowqual_Firm4_ids <- orig_Firm4_ids[-which(orig_Firm4_ids %in% checkm_output_highqual$accession)]
write.table(x = highqual_Firm4_ids, file = "highqual_genome_ids/Firm4_ids.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

orig_Firm5_ids <- read.table("all_genome_ids/Firm5_ids.txt", stringsAsFactors = FALSE)$V1
highqual_Firm5_ids <- orig_Firm5_ids[which(orig_Firm5_ids %in% checkm_output_highqual$accession)]
lowqual_Firm5_ids <- orig_Firm5_ids[-which(orig_Firm5_ids %in% checkm_output_highqual$accession)]
write.table(x = highqual_Firm5_ids, file = "highqual_genome_ids/Firm5_ids.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

orig_Gilliamella_ids <- read.table("all_genome_ids/Gilliamella_ids.txt", stringsAsFactors = FALSE)$V1
highqual_Gilliamella_ids <- orig_Gilliamella_ids[which(orig_Gilliamella_ids %in% checkm_output_highqual$accession)]
lowqual_Gilliamella_ids <- orig_Gilliamella_ids[-which(orig_Gilliamella_ids %in% checkm_output_highqual$accession)]
write.table(x = highqual_Gilliamella_ids, file = "highqual_genome_ids/Gilliamella_ids.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

orig_Snodgrasella_ids <- read.table("all_genome_ids/Snodgrasella_ids.txt", stringsAsFactors = FALSE)$V1
highqual_Snodgrasella_ids <- orig_Snodgrasella_ids[which(orig_Snodgrasella_ids %in% checkm_output_highqual$accession)]
lowqual_Snodgrasella_ids <- orig_Snodgrasella_ids[-which(orig_Snodgrasella_ids %in% checkm_output_highqual$accession)]
write.table(x = highqual_Snodgrasella_ids, file = "highqual_genome_ids/Snodgrasella_ids.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

