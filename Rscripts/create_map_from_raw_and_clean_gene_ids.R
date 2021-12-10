rm(list = ls(all.names = TRUE))

# Make mapfile from raw gene ids from prokka to ids with genome accessions in name.
# Also add phylotype name in front of gene ids that contain "refound".

Bifidobacterium_data <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Bifidobacterium/gene_data.csv",
                                   sep = ",", stringsAsFactors = FALSE, header = TRUE, quote = "", comment.char = "")

Bifidobacterium_data$gff_file <- gsub("GCF_", "GCF.", Bifidobacterium_data$gff_file)
Bifidobacterium_data$gff_file <- gsub("GCA_", "GCA.", Bifidobacterium_data$gff_file)
Bifidobacterium_data$gff_file <- gsub("_.*$", "", Bifidobacterium_data$gff_file)

Bifidobacterium_map <- data.frame(raw = Bifidobacterium_data$annotation_id,
                                  clean = paste(Bifidobacterium_data$gff_file, Bifidobacterium_data$annotation_id, sep = "_"))

Bifidobacterium_map[grep("refound", Bifidobacterium_map$raw), "raw"] <- paste("Bifidobacterium", Bifidobacterium_map[grep("refound", Bifidobacterium_map$raw), "raw"], sep = "_")



Firm4_data <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Firm4/gene_data.csv",
                                   sep = ",", stringsAsFactors = FALSE, header = TRUE, quote = "", comment.char = "")

Firm4_data$gff_file <- gsub("GCF_", "GCF.", Firm4_data$gff_file)
Firm4_data$gff_file <- gsub("GCA_", "GCA.", Firm4_data$gff_file)
Firm4_data$gff_file <- gsub("_.*$", "", Firm4_data$gff_file)

Firm4_map <- data.frame(raw = Firm4_data$annotation_id,
                                  clean = paste(Firm4_data$gff_file, Firm4_data$annotation_id, sep = "_"))

Firm4_map[grep("refound", Firm4_map$raw), "raw"] <- paste("Firm4", Firm4_map[grep("refound", Firm4_map$raw), "raw"], sep = "_")



Firm5_data <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Firm5/gene_data.csv",
                                   sep = ",", stringsAsFactors = FALSE, header = TRUE, quote = "", comment.char = "")

Firm5_data$gff_file <- gsub("GCF_", "GCF.", Firm5_data$gff_file)
Firm5_data$gff_file <- gsub("GCA_", "GCA.", Firm5_data$gff_file)
Firm5_data$gff_file <- gsub("_.*$", "", Firm5_data$gff_file)

Firm5_map <- data.frame(raw = Firm5_data$annotation_id,
                                  clean = paste(Firm5_data$gff_file, Firm5_data$annotation_id, sep = "_"))

Firm5_map[grep("refound", Firm5_map$raw), "raw"] <- paste("Firm5", Firm5_map[grep("refound", Firm5_map$raw), "raw"], sep = "_")



Gilliamella_data <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Gilliamella/gene_data.csv",
                                   sep = ",", stringsAsFactors = FALSE, header = TRUE, quote = "", comment.char = "")

Gilliamella_data$gff_file <- gsub("GCF_", "GCF.", Gilliamella_data$gff_file)
Gilliamella_data$gff_file <- gsub("GCA_", "GCA.", Gilliamella_data$gff_file)
Gilliamella_data$gff_file <- gsub("_.*$", "", Gilliamella_data$gff_file)

Gilliamella_map <- data.frame(raw = Gilliamella_data$annotation_id,
                                  clean = paste(Gilliamella_data$gff_file, Gilliamella_data$annotation_id, sep = "_"))

Gilliamella_map[grep("refound", Gilliamella_map$raw), "raw"] <- paste("Gilliamella", Gilliamella_map[grep("refound", Gilliamella_map$raw), "raw"], sep = "_")


Snodgrassella_data <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Snodgrassella/gene_data.csv",
                                   sep = ",", stringsAsFactors = FALSE, header = TRUE, quote = "", comment.char = "")

Snodgrassella_data$gff_file <- gsub("GCF_", "GCF.", Snodgrassella_data$gff_file)
Snodgrassella_data$gff_file <- gsub("GCA_", "GCA.", Snodgrassella_data$gff_file)
Snodgrassella_data$gff_file <- gsub("_.*$", "", Snodgrassella_data$gff_file)

Snodgrassella_map <- data.frame(raw = Snodgrassella_data$annotation_id,
                                  clean = paste(Snodgrassella_data$gff_file, Snodgrassella_data$annotation_id, sep = "_"))

Snodgrassella_map[grep("refound", Snodgrassella_map$raw), "raw"] <- paste("Snodgrassella", Snodgrassella_map[grep("refound", Snodgrassella_map$raw), "raw"], sep = "_")


combined_map <- do.call(rbind,
                        list(Bifidobacterium_map, Firm4_map, Firm5_map, Gilliamella_map, Snodgrassella_map))

write.table(x = combined_map,
            file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/map_raw_to_clean_gene_ids.tsv",
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)