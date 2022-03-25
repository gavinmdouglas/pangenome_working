### Get mean, sd, and median plasclass probability per ortholog (at least for species that are present).

rm(list = ls(all.names = TRUE))

panaroo_out <- list()

species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species_present.txt", stringsAsFactors = FALSE)$V1

all_orthologs <- c()

for (sp in species) {

  panaroo_out[[sp]] <- read.table(paste("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/", sp, "/gene_presence_absence.csv.gz", sep = ""),
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
  
  rownames(panaroo_out[[sp]]) <- paste(sp, panaroo_out[[sp]]$Gene, sep = "_")
  
  all_orthologs <- c(all_orthologs, rownames(panaroo_out[[sp]]))
}

plasclass_out <- read.table("/data1/gdouglas/projects/honey_bee/ref_genomes/plasclass_out/plasclass_output.tsv",
                            header = FALSE, sep = "\t", row.names = 1)
colnames(plasclass_out) <- "prob"


genes2contigs <- read.table("/data1/gdouglas/projects/honey_bee/ref_genomes/plasclass_out/contigs2genes.tsv",
                            header = FALSE, sep = "\t")
genes2contigs <- genes2contigs[which(genes2contigs$V2 != "."), ]
rownames(genes2contigs) <- genes2contigs$V2
genes2contigs <- genes2contigs[, -2, drop = FALSE]
colnames(genes2contigs) <- "contig"


ortholog_plasmid_prob <- data.frame(matrix(NA, nrow = length(all_orthologs), ncol = 5))
rownames(ortholog_plasmid_prob) <- all_orthologs
colnames(ortholog_plasmid_prob) <- c("gene", "max", "mean", "median", "sd")
ortholog_plasmid_prob$gene <- all_orthologs

for (sp in species) {
  
  sp_genes <- grep(sp, all_orthologs, value = TRUE)

  sp_num_columns <- ncol(panaroo_out[[sp]])
  
  for (g in sp_genes) {

    sp_gene_ids <- as.character(panaroo_out[[sp]][g, 4:sp_num_columns])

    # Sometimes gene ids are missing because they are "refound" (by panaroo), which are skipped for now.
    sp_gene_ids_matched <- sp_gene_ids[which(sp_gene_ids %in% rownames(genes2contigs))]
    
    num_matched_ids <- length(sp_gene_ids_matched)
    
    if (num_matched_ids == 0) { next }
    
    plasmid_prob <- plasclass_out[genes2contigs[sp_gene_ids_matched, "contig"], "prob"]
    
    prop_above_0.5 <- length(which(plasmid_prob > 0.5)) / length(plasmid_prob)
    
    ortholog_plasmid_prob[g, c("prop_above_0.5", "max", "mean", "median", "sd")] <- c(prop_above_0.5,
                                                                                      max(plasmid_prob),
                                                                                      mean(plasmid_prob),
                                                                                      median(plasmid_prob),
                                                                                      sd(plasmid_prob))
  
  }

}

write.table(x = ortholog_plasmid_prob,
            file = "/data1/gdouglas/projects/honey_bee/ref_genomes/plasclass_out/plasclass_output_ortholog_summary.tsv",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

