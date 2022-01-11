### Identify nonparalogs in the panaroo output

rm(list = ls(all.names = TRUE))

library(Biostrings)

# Remove certain gene ids not output in pangenome file.
Gilliamella_fasta_ids <- names(readDNAStringSet("honey_bee_pangenome/data/Ellegaard_panaroo_pangenome_mapping/clean_panaroo_pangenomes/Gilliamella_panaroo.fa"))

Gilliamella_fasta_ids <- gsub("Gilli_", "", Gilliamella_fasta_ids)

Gilliamella_panaroo <- read.table("honey_bee_pangenome/data/binning_output/panaroo_output_dastools_and_ref/Gilliamella_gene_presence_absence_roary-formatted.csv",
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")

Gilliamella_panaroo <- Gilliamella_panaroo[which(Gilliamella_panaroo$Gene %in% Gilliamella_fasta_ids), ]

Gilliamella_panaroo_clear_paralogs <- which(Gilliamella_panaroo$Avg.sequences.per.isolate > 1.1)

Gilliamella_panaroo <- Gilliamella_panaroo[-Gilliamella_panaroo_clear_paralogs, ]

##### Note that the "~~~" syntax just highlights cases where the gene was annotated differently in different prokka outputs - it doesn't necessarily mean they are paralogous!
##### Gilliamella_multifuncs <- grep("~~~", Gilliamella_panaroo$Gene)
##### Gilliamella_panaroo <- Gilliamella_panaroo[-Gilliamella_multifuncs, ]

Gilliamella_genes_to_keep <- paste("Gilli_", Gilliamella_panaroo$Gene, sep="")

write.table(file = "honey_bee_pangenome/data/Ellegaard_panaroo_pangenome_mapping/nonparalog_genes/Gilliamella_nonparalog.txt",
            x = Gilliamella_genes_to_keep,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


Snodgrassella_fasta_ids <- names(readDNAStringSet("honey_bee_pangenome/data/Ellegaard_panaroo_pangenome_mapping/clean_panaroo_pangenomes/Snodgrassella_panaroo.fa"))

Snodgrassella_fasta_ids <- gsub("Snod_", "", Snodgrassella_fasta_ids)

Snodgrassella_panaroo <- read.table("honey_bee_pangenome/data/binning_output/panaroo_output_dastools_and_ref/Snodgrassella_gene_presence_absence_roary-formatted.csv",
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")

Snodgrassella_panaroo <- Snodgrassella_panaroo[which(Snodgrassella_panaroo$Gene %in% Snodgrassella_fasta_ids), ]


Snodgrassella_panaroo_clear_paralogs <- which(Snodgrassella_panaroo$Avg.sequences.per.isolate > 1.1)

Snodgrassella_panaroo <- Snodgrassella_panaroo[-Snodgrassella_panaroo_clear_paralogs, ]

#####Snodgrassella_multifuncs <- grep("~~~", Snodgrassella_panaroo$Gene)
#####Snodgrassella_panaroo <- Snodgrassella_panaroo[-Snodgrassella_multifuncs, ]

Snodgrassella_genes_to_keep <- paste("Snod_", Snodgrassella_panaroo$Gene, sep="")

write.table(file = "honey_bee_pangenome/data/Ellegaard_panaroo_pangenome_mapping/nonparalog_genes/Snodgrassella_nonparalog.txt",
            x = Snodgrassella_genes_to_keep,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
