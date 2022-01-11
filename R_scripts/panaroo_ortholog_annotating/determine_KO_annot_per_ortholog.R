### Get annotations for overall orthologs.

rm(list = ls(all.names = TRUE))

library(stringr)

# Read in input
KO_output <- read.table(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/KEGG_output/all_microbiota_KO_calls.tsv.gz",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)

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


KO_calls <- list()
for (sp in species) {
  
  print(sp)
  
  for (ortholog in rownames(panaroo_out[[sp]])) {
    ortholog_KO_calls <- c()
    num_CDS <- 0
    for (calls in panaroo_out[[sp]][ortholog, 4:ncol(panaroo_out[[sp]])]) {
         for (call_id in str_split(calls, ";")[[1]]) {
          if ((call_id == "") | is.na(call_id)) { next }
            num_CDS <- num_CDS + 1
            
            # Remove suffixes indicating short lengths or premature stop codons.
            call_id <- gsub("_len$", "", call_id)
            call_id <- gsub("_stop$", "", call_id)
            
            if (call_id %in% rownames(KO_output)) {
              ortholog_KO_calls <- c(ortholog_KO_calls, str_split(KO_output[call_id, "KO"], ",")[[1]])
            }
       }
    }
    
    if (length(ortholog_KO_calls) > 0) {
      ortholog_KO_calls_prop <- table(ortholog_KO_calls) / num_CDS
      
      if (length(which(ortholog_KO_calls_prop >= 0.9)) > 0) {
        KO_calls[[ortholog]] <- paste(sort(names(ortholog_KO_calls_prop[which(ortholog_KO_calls_prop >= 0.9)])), collapse = ",")
      }
    }
  }
}

KO_calls_clean <- do.call(c, KO_calls)

KO_calls_out <- data.frame(gene=names(KO_calls_clean), KO=KO_calls_clean)

write.table(x = KO_calls_out,
            file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/KEGG_output/all_microbiota_KO_by_ortholog.tsv",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


KO_call_tally_tab <- data.frame(matrix(NA, nrow = length(species), ncol = 2))
colnames(KO_call_tally_tab) <- c("count", "percent")
rownames(KO_call_tally_tab) <- species

for (sp in rownames(KO_call_tally_tab)) {
  KO_call_tally_tab[sp, ] <- c(length(grep(sp, KO_calls_out$gene)),
                                        (length(grep(sp, KO_calls_out$gene)) / nrow(panaroo_out[[sp]])) * 100)
}

write.table(x = KO_call_tally_tab,
            file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/KEGG_output/all_microbiota_KO_by_ortholog_species_tallies.tsv",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")
