rm(list = ls(all.names = TRUE))

library("discover")
library("plyr")

source("projects/honey_bee/scripts/pangenome_working/Rscripts/functions.R")


### Gilliamella ###

Gilliamella_path_to_panaroo <- "/home/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Gilliamella_gene_presence_absence_roary-formatted.csv.gz"
Gilliamella_path_to_depth_Ellegaard.2019 <- "/home/gdouglas/projects/honey_bee/Ellegaard/pangenome_competitive_mapping/bedtools_coverage_out/panaroo_best_hit/Gilliamella_mean_depth_per_site/"
Gilliamella_path_to_depth_Ellegaard.2020 <- "/data1/gdouglas/projects/honey_bee/Ellegaard.2020/panaroo_bedtools_coverage/Gilliamella_mean_depth_per_site/"

Snodgrassella_path_to_panaroo <- "/home/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Snodgrassella_gene_presence_absence_roary-formatted.csv.gz"
Snodgrassella_path_to_depth_Ellegaard.2019  <- "/home/gdouglas/projects/honey_bee/Ellegaard/pangenome_competitive_mapping/bedtools_coverage_out/panaroo_best_hit/Snodgrassella_mean_depth_per_site/"
Snodgrassella_path_to_depth_Ellegaard.2020 <- "/data1/gdouglas/projects/honey_bee/Ellegaard.2020/panaroo_bedtools_coverage/Snodgrassella_mean_depth_per_site/"

Gilliamella_panaroo_out <- read.table(Gilliamella_path_to_panaroo,
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilli", rownames(Gilliamella_panaroo_out), sep = "_")


Gilliamella_basic_depth_Ellegaard.2019 <- read_in_depth_files(in_path = Gilliamella_path_to_depth_Ellegaard.2019,
                                                                  pattern = ".mean.bedGraph.gz",
                                                                  roary_formatted_pangenome = Gilliamella_panaroo_out)

Gilliamella_basic_depth_Ellegaard.2020 <- read_in_depth_files(in_path = Gilliamella_path_to_depth_Ellegaard.2020,
                                                                  pattern = ".mean.bedGraph.gz",
                                                                  roary_formatted_pangenome = Gilliamella_panaroo_out)

Gilliamella_depth_by_sample_clean <- list()

Gilliamella_Ellegaard.2019_samples <- names(Gilliamella_basic_depth_Ellegaard.2019$depth_by_sample)
Gilliamella_Ellegaard.2020_samples <- names(Gilliamella_basic_depth_Ellegaard.2020$depth_by_sample)

for (Gilliamella_Ellegaard.2019_sample in Gilliamella_Ellegaard.2019_samples) {
  Gilliamella_depth_by_sample_clean[[Gilliamella_Ellegaard.2019_sample]] <- Gilliamella_basic_depth_Ellegaard.2019$depth_by_sample[[Gilliamella_Ellegaard.2019_sample]][, c("gene", "V4"), drop = FALSE]
  colnames(Gilliamella_depth_by_sample_clean[[Gilliamella_Ellegaard.2019_sample]]) <- c("gene", Gilliamella_Ellegaard.2019_sample)
}

for (Gilliamella_Ellegaard.2020_sample in Gilliamella_Ellegaard.2020_samples) {
  Gilliamella_depth_by_sample_clean[[Gilliamella_Ellegaard.2020_sample]] <- Gilliamella_basic_depth_Ellegaard.2020$depth_by_sample[[Gilliamella_Ellegaard.2020_sample]][, c("gene", "V4"), drop = FALSE]
  colnames(Gilliamella_depth_by_sample_clean[[Gilliamella_Ellegaard.2020_sample]]) <- c("gene", Gilliamella_Ellegaard.2020_sample)
}


Gilliamella_depth <- join_all(Gilliamella_depth_by_sample_clean, by="gene", type='left')
rownames(Gilliamella_depth) <- Gilliamella_depth$gene
Gilliamella_depth <- Gilliamella_depth[, -which(colnames(Gilliamella_depth) == "gene")]

saveRDS(object = Gilliamella_depth,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_mean_depth_per_site.rds")


### Snodgrassella ###

Snodgrassella_panaroo_out <- read.table(Snodgrassella_path_to_panaroo,
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snod", rownames(Snodgrassella_panaroo_out), sep = "_")


Snodgrassella_basic_depth_Ellegaard.2019 <- read_in_depth_files(in_path = Snodgrassella_path_to_depth_Ellegaard.2019,
                                                              pattern = ".mean.bedGraph.gz",
                                                              roary_formatted_pangenome = Snodgrassella_panaroo_out)

Snodgrassella_basic_depth_Ellegaard.2020 <- read_in_depth_files(in_path = Snodgrassella_path_to_depth_Ellegaard.2020,
                                                              pattern = ".mean.bedGraph.gz",
                                                              roary_formatted_pangenome = Snodgrassella_panaroo_out)

Snodgrassella_depth_by_sample_clean <- list()

Snodgrassella_Ellegaard.2019_samples <- names(Snodgrassella_basic_depth_Ellegaard.2019$depth_by_sample)
Snodgrassella_Ellegaard.2020_samples <- names(Snodgrassella_basic_depth_Ellegaard.2020$depth_by_sample)

for (Snodgrassella_Ellegaard.2019_sample in Snodgrassella_Ellegaard.2019_samples) {
  Snodgrassella_depth_by_sample_clean[[Snodgrassella_Ellegaard.2019_sample]] <- Snodgrassella_basic_depth_Ellegaard.2019$depth_by_sample[[Snodgrassella_Ellegaard.2019_sample]][, c("gene", "V4"), drop = FALSE]
  colnames(Snodgrassella_depth_by_sample_clean[[Snodgrassella_Ellegaard.2019_sample]]) <- c("gene", Snodgrassella_Ellegaard.2019_sample)
}

for (Snodgrassella_Ellegaard.2020_sample in Snodgrassella_Ellegaard.2020_samples) {
  Snodgrassella_depth_by_sample_clean[[Snodgrassella_Ellegaard.2020_sample]] <- Snodgrassella_basic_depth_Ellegaard.2020$depth_by_sample[[Snodgrassella_Ellegaard.2020_sample]][, c("gene", "V4"), drop = FALSE]
  colnames(Snodgrassella_depth_by_sample_clean[[Snodgrassella_Ellegaard.2020_sample]]) <- c("gene", Snodgrassella_Ellegaard.2020_sample)
}


Snodgrassella_depth <- join_all(Snodgrassella_depth_by_sample_clean, by="gene", type='left')
rownames(Snodgrassella_depth) <- Snodgrassella_depth$gene
Snodgrassella_depth <- Snodgrassella_depth[, -which(colnames(Snodgrassella_depth) == "gene")]

saveRDS(object = Snodgrassella_depth,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Snodgrassella_mean_depth_per_site.rds")

