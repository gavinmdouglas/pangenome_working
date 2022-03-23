### Get RDS of list of categories to genes in each set, which is convenient for running functional enrichment.

rm(list = ls(all.names = TRUE))

ortholog_to_category <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/all_microbiota_COG.category_by_ortholog.tsv",
                                   header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

all_categories <- read.table("/data1/gdouglas/db/COG_definitions/COG_category_descrip.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

category_to_genes <- list()

for (category in all_categories) {
  
  category_genes <- rownames(ortholog_to_category)[grep(category, ortholog_to_category$V2)]
  
  if (length(category_genes) > 0) {
    category_to_genes[[category]] <- rownames(ortholog_to_category)[grep(category, ortholog_to_category$V2)]
  }
}

category_to_genes[["any_category"]] <- rownames(ortholog_to_category)

saveRDS(object = category_to_genes, file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/COG_category_to_ortholog.rds")