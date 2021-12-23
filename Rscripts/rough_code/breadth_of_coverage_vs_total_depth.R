rm(list = ls(all.names = TRUE))

source("honey_bee_pangenome/Rscripts/functions.R")

# Read in reference genome info.
Gilliamella_panaroo_out <- read.table("honey_bee_pangenome/data/binning_output/panaroo_output_dastools_and_ref/Gilliamella_gene_presence_absence_roary-formatted.csv",
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilli", rownames(Gilliamella_panaroo_out), sep = "_")

Gilliamella_basic_breadth <- read_in_breadth_files(in_path = "honey_bee_pangenome/data/Ellegaard_panaroo_pangenome_mapping/bowtie2_mapped_gene_coverage/best_hit_simple_comp_map/Gilliamella_coverage_breadth/",
                                                   pattern = ".breadth.bedGraph.gz",
                                                   roary_formatted_pangenome = Gilliamella_panaroo_out)

Gilliamella_basic_panaroo_core <- rownames(Gilliamella_panaroo_out)[which(Gilliamella_panaroo_out$No..isolates >= 99)]

Gilliamella_mean_core_breadth <- c()
Gilliamella_median_core_breadth <- c()

for (SRR in names(Gilliamella_basic_breadth$breadth_by_sample)) {

  Gilliamella_mean_core_breadth <- c(Gilliamella_mean_core_breadth,
                                     mean(Gilliamella_basic_breadth$breadth_by_sample[[SRR]][which(Gilliamella_basic_breadth$breadth_by_sample[[SRR]]$gene %in% Gilliamella_basic_panaroo_core), "V7"]))
 
   Gilliamella_median_core_breadth <- c(Gilliamella_median_core_breadth,
                                     median(Gilliamella_basic_breadth$breadth_by_sample[[SRR]][which(Gilliamella_basic_breadth$breadth_by_sample[[SRR]]$gene %in% Gilliamella_basic_panaroo_core), "V7"]))
  
}

names(Gilliamella_mean_core_breadth) <- names(Gilliamella_basic_breadth$breadth_by_sample)
names(Gilliamella_median_core_breadth) <- names(Gilliamella_basic_breadth$breadth_by_sample)

Gilliamella_mapped_read_counts_df <- read.table("honey_bee_pangenome/data/tests/test_coverage_breadth_no_competitive_mapping/panaroo_best_hit_mapped_read_counts/Gilliamella_nonparalog_read_counts.txt",
                                             header = FALSE, sep =" ", stringsAsFactors = FALSE)

Gilliamella_mapped_read_counts <- Gilliamella_mapped_read_counts_df[, 2]
names(Gilliamella_mapped_read_counts) <- Gilliamella_mapped_read_counts_df[, 1]

Gilliamella_mapped_read_counts <- Gilliamella_mapped_read_counts[names(Gilliamella_mean_core_breadth)]

plot(Gilliamella_mapped_read_counts, Gilliamella_median_core_breadth)

# Read in reference genome info.
Snodgrassella_panaroo_out <- read.table("honey_bee_pangenome/data/binning_output/panaroo_output_dastools_and_ref/Snodgrassella_gene_presence_absence_roary-formatted.csv",
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snod", rownames(Snodgrassella_panaroo_out), sep = "_")

Snodgrassella_basic_breadth <- read_in_breadth_files(in_path = "honey_bee_pangenome/data/Ellegaard_panaroo_pangenome_mapping/bowtie2_mapped_gene_coverage/best_hit_simple_comp_map/Snodgrassella_coverage_breadth/",
                                                     pattern = ".breadth.bedGraph.gz",
                                                     roary_formatted_pangenome = Snodgrassella_panaroo_out)

Snodgrassella_basic_panaroo_core <- rownames(Snodgrassella_panaroo_out)[which(Snodgrassella_panaroo_out$No..isolates >= 53)]


Snodgrassella_mean_core_breadth <- c()
Snodgrassella_median_core_breadth <- c()

for (SRR in names(Snodgrassella_basic_breadth$breadth_by_sample)) {
  
  Snodgrassella_mean_core_breadth <- c(Snodgrassella_mean_core_breadth,
                                     mean(Snodgrassella_basic_breadth$breadth_by_sample[[SRR]][which(Snodgrassella_basic_breadth$breadth_by_sample[[SRR]]$gene %in% Snodgrassella_basic_panaroo_core), "V7"]))
  
  Snodgrassella_median_core_breadth <- c(Snodgrassella_median_core_breadth,
                                       median(Snodgrassella_basic_breadth$breadth_by_sample[[SRR]][which(Snodgrassella_basic_breadth$breadth_by_sample[[SRR]]$gene %in% Snodgrassella_basic_panaroo_core), "V7"]))
  
}

names(Snodgrassella_mean_core_breadth) <- names(Snodgrassella_basic_breadth$breadth_by_sample)
names(Snodgrassella_median_core_breadth) <- names(Snodgrassella_basic_breadth$breadth_by_sample)

Snodgrassella_mapped_read_counts_df <- read.table("honey_bee_pangenome/data/tests/test_coverage_breadth_no_competitive_mapping/panaroo_best_hit_mapped_read_counts/Snodgrassella_nonparalog_read_counts.txt",
                                                header = FALSE, sep =" ", stringsAsFactors = FALSE)

Snodgrassella_mapped_read_counts <- Snodgrassella_mapped_read_counts_df[, 2]
names(Snodgrassella_mapped_read_counts) <- Snodgrassella_mapped_read_counts_df[, 1]

Snodgrassella_mapped_read_counts <- Snodgrassella_mapped_read_counts[names(Snodgrassella_mean_core_breadth)]

plot(Snodgrassella_mapped_read_counts, Snodgrassella_median_core_breadth)
