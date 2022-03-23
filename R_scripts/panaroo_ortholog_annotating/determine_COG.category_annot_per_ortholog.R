### Get annotations for overall orthologs.
### Note that the COG category classification in the eggNOG output is based on an outdated version that doesn't include mobile elements.
### So I re-ran the mapping based on cases where the broadest OG was a COG.

rm(list = ls(all.names = TRUE))

library(stringr)

# Read in input
eggNOG_output <- read.table(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/all_microbiota_eggNOG.tsv",
                            header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, quote = "", comment.char = "")

COG2category <- readRDS("/data1/gdouglas/db/COG_definitions/cog-20.to_category_collapse.rds")

COG_output <- eggNOG_output[grep("^COG", eggNOG_output$broadest_OG), ]

# Noted that a few of the COGs were missing from the COG to category mapping, so I added the mapping manually.
# These were based on the definitions on NCBI for the gene family
# E.g., https://www.ncbi.nlm.nih.gov/Structure/cdd/COG5061
unique(COG_output$broadest_OG[which(! COG_output$broadest_OG %in% rownames(COG2category))])
COG2category["COG3512", ] <- c("COG3512", "V")
COG2category["COG5290", ] <- c("COG5290", "K")
COG2category["COG5113", ] <- c("COG5113", "O")
COG2category["COG5061", ] <- c("COG5061", "O,U")

COG_output$COG_category <- COG2category[COG_output$broadest_OG, "category"]

# # Need to add commands between cases where genes are assigned to multiple COG categories.
# NOTE: These commands would only be needed if you were to use the original eggNOG COG category column, which does not delimit by commas
# COG_category_nchar <- sapply(COG_output$COG_category, nchar)
# 
# for (i in which(COG_category_nchar != 1)) {
#   
#   new_COG_category <- paste(sort(strsplit(COG_output[i, "COG_category"], split = "")[[1]]), collapse = ",")
#   
#   COG_output[i, "COG_category"] <- new_COG_category
# }

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


COG_calls <- list()
for (sp in species) {
  
  print(sp)
  
  for (ortholog in rownames(panaroo_out[[sp]])) {
    ortholog_COG_calls <- c()
    num_CDS <- 0
    for (calls in panaroo_out[[sp]][ortholog, 4:ncol(panaroo_out[[sp]])]) {
         for (call_id in str_split(calls, ";")[[1]]) {
          if ((call_id == "") | is.na(call_id)) { next }
            num_CDS <- num_CDS + 1
            
            # Remove suffixes indicating short lengths or premature stop codons.
            call_id <- gsub("_len$", "", call_id)
            call_id <- gsub("_stop$", "", call_id)
            
            if (call_id %in% rownames(COG_output)) {
              ortholog_COG_calls <- c(ortholog_COG_calls, str_split(COG_output[call_id, "COG_category"], ",")[[1]])
            }
       }
    }
    
    if (length(ortholog_COG_calls) > 0) {
      ortholog_COG_calls_prop <- table(ortholog_COG_calls) / num_CDS
      
      if (length(which(ortholog_COG_calls_prop >= 0.9)) > 0) {
        COG_calls[[ortholog]] <- paste(sort(names(ortholog_COG_calls_prop[which(ortholog_COG_calls_prop >= 0.9)])), collapse = ",")
      }
    }
  }
}

COG_calls_clean <- do.call(c, COG_calls)

COG_calls_out <- data.frame(gene=names(COG_calls_clean), COG=COG_calls_clean)

write.table(x = COG_calls_out,
            file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/all_microbiota_COG.category_by_ortholog.tsv",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


COG_call_tally_tab <- data.frame(matrix(NA, nrow = length(species), ncol = 2))
colnames(COG_call_tally_tab) <- c("count", "percent")
rownames(COG_call_tally_tab) <- species

for (sp in rownames(COG_call_tally_tab)) {
  COG_call_tally_tab[sp, ] <- c(length(grep(sp, COG_calls_out$gene)),
                                        (length(grep(sp, COG_calls_out$gene)) / nrow(panaroo_out[[sp]])) * 100)
}

write.table(x = COG_call_tally_tab,
            file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/all_microbiota_COG.category_by_ortholog_species_tallies.tsv",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")
