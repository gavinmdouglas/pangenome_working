### Merge eggNOG-mapper output tables into single table.

rm(list = ls(all.names = TRUE))

eggNOG <- list()
eggNOG[["part1"]] <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/raw/all_microbiota_proteins.part-1_eggNOG_mapper_out.txt.gz",
                                header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, comment.char = "", quote = "")
eggNOG[["part2"]] <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/raw/all_microbiota_proteins.part-2_eggNOG_mapper_out.txt.gz",
                                header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, comment.char = "", quote = "")
eggNOG[["part3"]] <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/raw/all_microbiota_proteins.part-3_eggNOG_mapper_out.txt.gz",
                                header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, comment.char = "", quote = "")
eggNOG[["part4"]] <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/raw/all_microbiota_proteins.part-4_eggNOG_mapper_out.txt.gz",
                                header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, comment.char = "", quote = "")
eggNOG[["part5"]] <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/raw/all_microbiota_proteins.part-5_eggNOG_mapper_out.txt.gz",
                                header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, comment.char = "", quote = "")
eggNOG_all_hits <- do.call(what = rbind, args = eggNOG)

rownames(eggNOG_all_hits) <- gsub("^part\\d\\.", "", rownames(eggNOG_all_hits))

eggNOG_all_hits[eggNOG_all_hits == "-"] <- NA

# Remove all KEGG info except for KOs (in case sanity checks are needed later)
eggNOG_all_hits <- eggNOG_all_hits[, -which(colnames(eggNOG_all_hits) %in% c("KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC"))]

eggNOG_all_hits$KEGG_ko <- gsub("^ko:", "", eggNOG_all_hits$KEGG_ko)


# Clean up eggNOG id column to get broadest OG.
eggNOG_all_hits$broadest_OG <- gsub("@.*$", "", eggNOG_all_hits$eggNOG_OGs)
eggNOG_all_hits$broadest_OG <- gsub("\"", "", eggNOG_all_hits$broadest_OG)


# Percent of CDS sequences annotated with eggNOG.
num_CDS <- 439817
(nrow(eggNOG_all_hits) / num_CDS) * 100
(length(which(! is.na(eggNOG_all_hits$COG_category))) / num_CDS) * 100
(length(which(! is.na(eggNOG_all_hits$CAZy))) / num_CDS) * 100

write.table(x = eggNOG_all_hits,
            file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/all_microbiota_eggNOG.tsv",
            row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")

