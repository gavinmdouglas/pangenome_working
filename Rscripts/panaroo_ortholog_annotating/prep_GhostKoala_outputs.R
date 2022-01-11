### Post-process raw GhostKoala output

# First had to run this command to only keep rows where a KO was called.
# awk '{ if ( NF == 2 ) print $0  }' user_ko.txt > all_microbiota_KO_calls_long.tsv

rm(list = ls(all.names = TRUE))

KO_raw_out <- read.table(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/KEGG_output/all_microbiota_KO_calls_long.tsv.gz",
                         header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

uniq_genes <- sort(unique(KO_raw_out$V1))


KO_clean_out <- data.frame(matrix(NA, nrow = length(uniq_genes), ncol = 2))
colnames(KO_clean_out) <- c("CDS", "KO")
KO_clean_out$CDS <- uniq_genes
rownames(KO_clean_out) <- uniq_genes

for (gene in uniq_genes) {
  KO_hits <- KO_raw_out[which(KO_raw_out$V1 == gene), "V2"]
  
  if (length(KO_hits) == 1) {
    KO_clean_out[gene, "KO"] <- KO_hits
  } else {
    KO_clean_out[gene, "KO"] <- paste(sort(KO_hits), collapse = ",")
  }
}



write.table(x = KO_clean_out,
            file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/KEGG_output/all_microbiota_KO_calls.tsv",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

