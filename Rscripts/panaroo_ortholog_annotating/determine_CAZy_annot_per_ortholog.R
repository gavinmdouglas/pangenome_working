### Get annotations for overall orthologs.

rm(list = ls(all.names = TRUE))

library(stringr)

# Read in input
eggNOG_output <- read.table(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/all_microbiota_eggNOG.tsv",
                            header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, quote = "", comment.char = "")

CAZy_output <- eggNOG_output[which(!is.na(eggNOG_output$CAZy)), ]

panaroo_out <- list()

species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species.txt", stringsAsFactors = FALSE)$V1

for (sp in species) {
  
  if (sp == "Apilactobacillus_apinorum") { next } # Only 1 genome so treat differently.
  
  panaroo_out[[sp]] <- read.table(paste("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/", sp, "/gene_presence_absence.csv.gz", sep = ""),
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
  
  rownames(panaroo_out[[sp]]) <- paste(sp, panaroo_out[[sp]]$Gene, sep = "_")
  
}

# Create dummy Apilactobacillus_apinorum table.
Apilactobacillus_apinorum_ids <- read.table("/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes_prokka/Apilactobacillus_apinorum/GCA_001281175.1/GCA_001281175.1.ids.txt",
                                            header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1
panaroo_out[["Apilactobacillus_apinorum"]] <- data.frame(matrix(NA, nrow = length(Apilactobacillus_apinorum_ids), ncol = 4))
colnames(panaroo_out[["Apilactobacillus_apinorum"]]) <- c("Gene", "Non.unique.Gene.name", "Annotation", "GCA_001281175.1")
rownames(panaroo_out[["Apilactobacillus_apinorum"]]) <- paste("Apilactobacillus_apinorum", Apilactobacillus_apinorum_ids, sep = "_")
panaroo_out[["Apilactobacillus_apinorum"]]$GCA_001281175.1 <- Apilactobacillus_apinorum_ids


CAZy_calls <- list()
for (sp in species) {
  
  print(sp)
  
  for (ortholog in rownames(panaroo_out[[sp]])) {
    ortholog_CAZy_calls <- c()
    num_CDS <- 0
    for (calls in panaroo_out[[sp]][ortholog, 4:ncol(panaroo_out[[sp]])]) {
         for (call_id in str_split(calls, ";")[[1]]) {
          if ((call_id == "") | is.na(call_id)) { next }
            num_CDS <- num_CDS + 1
            
            # Remove suffixes indicating short lengths or premature stop codons.
            call_id <- gsub("_len$", "", call_id)
            call_id <- gsub("_stop$", "", call_id)
            
            if (call_id %in% rownames(CAZy_output)) {
              ortholog_CAZy_calls <- c(ortholog_CAZy_calls, str_split(CAZy_output[call_id, "CAZy"], ",")[[1]])
            }
       }
    }
    
    if (length(ortholog_CAZy_calls) > 0) {
      ortholog_CAZy_calls_prop <- table(ortholog_CAZy_calls) / num_CDS
      
      if (length(which(ortholog_CAZy_calls_prop >= 0.9)) > 0) {
        CAZy_calls[[ortholog]] <- paste(sort(names(ortholog_CAZy_calls_prop[which(ortholog_CAZy_calls_prop >= 0.9)])), collapse = ",")
      }
    }
  }
}

CAZy_calls_clean <- do.call(c, CAZy_calls)

CAZy_calls_out <- data.frame(gene=names(CAZy_calls_clean), CAZy=CAZy_calls_clean)

write.table(x = CAZy_calls_out,
            file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/all_microbiota_CAZy_by_ortholog.tsv",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


CAZy_call_tally_tab <- data.frame(matrix(NA, nrow = length(species), ncol = 2))
colnames(CAZy_call_tally_tab) <- c("count", "percent")
rownames(CAZy_call_tally_tab) <- species

for (sp in rownames(CAZy_call_tally_tab)) {
  CAZy_call_tally_tab[sp, ] <- c(length(grep(sp, CAZy_calls_out$gene)),
                                        (length(grep(sp, CAZy_calls_out$gene)) / nrow(panaroo_out[[sp]])) * 100)
}

write.table(x = CAZy_call_tally_tab,
            file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/all_microbiota_CAZy_by_ortholog_species_tallies.tsv",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")
