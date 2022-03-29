### Prep presence absence matrices of both Snodgrassella and Gilliamella (for both Ellegaard 2020 and 2019) including noncore genes and not.

rm(list = ls(all.names = TRUE))

library("discover")
library("plyr")

source("projects/honey_bee/scripts/pangenome_working/Rscripts/functions.R")

Gilliamella_path_to_panaroo <- "/home/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Gilliamella_gene_presence_absence_roary-formatted.csv.gz"
Gilliamella_path_to_basic_breadth_Ellegaard.2019 <- "/home/gdouglas/projects/honey_bee/Ellegaard/pangenome_competitive_mapping/bedtools_coverage_out/panaroo_best_hit/Gilliamella_coverage_breadth/"
Gilliamella_path_to_basic_breadth_Ellegaard.2020 <- "/data1/gdouglas/projects/honey_bee/Ellegaard.2020/panaroo_bedtools_coverage/Gilliamella_coverage_breadth/"

Snodgrassella_path_to_panaroo <- "/home/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Snodgrassella_gene_presence_absence_roary-formatted.csv.gz"
Snodgrassella_path_to_basic_breadth_Ellegaard.2019  <- "/home/gdouglas/projects/honey_bee/Ellegaard/pangenome_competitive_mapping/bedtools_coverage_out/panaroo_best_hit/Snodgrassella_coverage_breadth/"
Snodgrassella_path_to_basic_breadth_Ellegaard.2020 <- "/data1/gdouglas/projects/honey_bee/Ellegaard.2020/panaroo_bedtools_coverage/Snodgrassella_coverage_breadth/"

Gilliamella_panaroo_out <- read.table(Gilliamella_path_to_panaroo,
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilli", rownames(Gilliamella_panaroo_out), sep = "_")


Gilliamella_basic_breadth_Ellegaard.2019 <- read_in_breadth_files(in_path = Gilliamella_path_to_basic_breadth_Ellegaard.2019,
                                                                  pattern = ".breadth.bedGraph.gz",
                                                                  roary_formatted_pangenome = Gilliamella_panaroo_out)

Gilliamella_basic_breadth_Ellegaard.2020 <- read_in_breadth_files(in_path = Gilliamella_path_to_basic_breadth_Ellegaard.2020,
                                                                  pattern = ".breadth.bedGraph.gz",
                                                                  roary_formatted_pangenome = Gilliamella_panaroo_out)

Gilliamella_breadth_by_sample_clean <- list()

Gilliamella_Ellegaard.2019_samples <- names(Gilliamella_basic_breadth_Ellegaard.2019$breadth_by_sample)
Gilliamella_Ellegaard.2020_samples <- names(Gilliamella_basic_breadth_Ellegaard.2020$breadth_by_sample)

for (Gilliamella_Ellegaard.2019_sample in Gilliamella_Ellegaard.2019_samples) {
  Gilliamella_breadth_by_sample_clean[[Gilliamella_Ellegaard.2019_sample]] <- Gilliamella_basic_breadth_Ellegaard.2019$breadth_by_sample[[Gilliamella_Ellegaard.2019_sample]][, c("gene", "V7"), drop = FALSE]
  colnames(Gilliamella_breadth_by_sample_clean[[Gilliamella_Ellegaard.2019_sample]]) <- c("gene", Gilliamella_Ellegaard.2019_sample)
}

for (Gilliamella_Ellegaard.2020_sample in Gilliamella_Ellegaard.2020_samples) {
  Gilliamella_breadth_by_sample_clean[[Gilliamella_Ellegaard.2020_sample]] <- Gilliamella_basic_breadth_Ellegaard.2020$breadth_by_sample[[Gilliamella_Ellegaard.2020_sample]][, c("gene", "V7"), drop = FALSE]
  colnames(Gilliamella_breadth_by_sample_clean[[Gilliamella_Ellegaard.2020_sample]]) <- c("gene", Gilliamella_Ellegaard.2020_sample)
}


Gilliamella_present <- join_all(Gilliamella_breadth_by_sample_clean, by="gene", type='left')
rownames(Gilliamella_present) <- Gilliamella_present$gene
Gilliamella_present <- Gilliamella_present[, -which(colnames(Gilliamella_present) == "gene")]
Gilliamella_present[Gilliamella_present < presence_cutoff] <- 0
Gilliamella_present[Gilliamella_present >= presence_cutoff] <- 1

Gilliamella_noncore_genes <- rownames(Gilliamella_panaroo_out)[which(Gilliamella_panaroo_out$No..isolates < 0.95 * max(Gilliamella_panaroo_out$No..isolates))]

Gilliamella_present_noncore <- Gilliamella_present[which(rownames(Gilliamella_present) %in% Gilliamella_noncore_genes), ]


Snodgrassella_panaroo_out <- read.table(Snodgrassella_path_to_panaroo,
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snod", rownames(Snodgrassella_panaroo_out), sep = "_")


Snodgrassella_basic_breadth_Ellegaard.2019 <- read_in_breadth_files(in_path = Snodgrassella_path_to_basic_breadth_Ellegaard.2019,
                                                                  pattern = ".breadth.bedGraph.gz",
                                                                  roary_formatted_pangenome = Snodgrassella_panaroo_out)

Snodgrassella_basic_breadth_Ellegaard.2020 <- read_in_breadth_files(in_path = Snodgrassella_path_to_basic_breadth_Ellegaard.2020,
                                                                  pattern = ".breadth.bedGraph.gz",
                                                                  roary_formatted_pangenome = Snodgrassella_panaroo_out)

Snodgrassella_breadth_by_sample_clean <- list()

Snodgrassella_Ellegaard.2019_samples <- names(Snodgrassella_basic_breadth_Ellegaard.2019$breadth_by_sample)
Snodgrassella_Ellegaard.2020_samples <- names(Snodgrassella_basic_breadth_Ellegaard.2020$breadth_by_sample)

for (Snodgrassella_Ellegaard.2019_sample in Snodgrassella_Ellegaard.2019_samples) {
  Snodgrassella_breadth_by_sample_clean[[Snodgrassella_Ellegaard.2019_sample]] <- Snodgrassella_basic_breadth_Ellegaard.2019$breadth_by_sample[[Snodgrassella_Ellegaard.2019_sample]][, c("gene", "V7"), drop = FALSE]
  colnames(Snodgrassella_breadth_by_sample_clean[[Snodgrassella_Ellegaard.2019_sample]]) <- c("gene", Snodgrassella_Ellegaard.2019_sample)
}

for (Snodgrassella_Ellegaard.2020_sample in Snodgrassella_Ellegaard.2020_samples) {
  Snodgrassella_breadth_by_sample_clean[[Snodgrassella_Ellegaard.2020_sample]] <- Snodgrassella_basic_breadth_Ellegaard.2020$breadth_by_sample[[Snodgrassella_Ellegaard.2020_sample]][, c("gene", "V7"), drop = FALSE]
  colnames(Snodgrassella_breadth_by_sample_clean[[Snodgrassella_Ellegaard.2020_sample]]) <- c("gene", Snodgrassella_Ellegaard.2020_sample)
}


Snodgrassella_present <- join_all(Snodgrassella_breadth_by_sample_clean, by="gene", type='left')
rownames(Snodgrassella_present) <- Snodgrassella_present$gene
Snodgrassella_present <- Snodgrassella_present[, -which(colnames(Snodgrassella_present) == "gene")]
Snodgrassella_present[Snodgrassella_present < presence_cutoff] <- 0
Snodgrassella_present[Snodgrassella_present >= presence_cutoff] <- 1

Snodgrassella_noncore_genes <- rownames(Snodgrassella_panaroo_out)[which(Snodgrassella_panaroo_out$No..isolates < 0.95 * max(Snodgrassella_panaroo_out$No..isolates))]

Snodgrassella_present_noncore <- Snodgrassella_present[which(rownames(Snodgrassella_present) %in% Snodgrassella_noncore_genes), ]


Gilliamella_and_Snodgrassella_present <- rbind(Gilliamella_present, Snodgrassella_present)
Gilliamella_and_Snodgrassella_present_noncore <- rbind(Gilliamella_present_noncore, Snodgrassella_present_noncore)

saveRDS(object = Gilliamella_and_Snodgrassella_present,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_and_Snodgrassella_present.rds")

saveRDS(object = Gilliamella_and_Snodgrassella_present_noncore,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_and_Snodgrassella_present_noncore.rds")

