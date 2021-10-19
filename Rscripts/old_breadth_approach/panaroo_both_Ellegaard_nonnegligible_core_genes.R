### Identify core genes that are consistently mapped to across samples in both datasets.

rm(list = ls(all.names = TRUE))

source("/home/gdouglas/projects/honey_bee/scripts/pangenome_working/Rscripts/functions.R")

Gilliamella_path_to_panaroo <- "/home/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Gilliamella_gene_presence_absence_roary-formatted.csv.gz"
Gilliamella_path_to_basic_breadth_combined <- "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/per_sample_tables/Gilliamella_coverage_breadth"

Snodgrassella_path_to_panaroo <- "/home/gdouglas/projects/honey_bee/panaroo_pangenomes/roary_formatted_files/Snodgrassella_gene_presence_absence_roary-formatted.csv.gz"
Snodgrassella_path_to_basic_breadth_combined <- "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/per_sample_tables/Snodgrassella_coverage_breadth"



### Gilliamella

# Read in reference genome info.
Gilliamella_panaroo_out <- read.table(Gilliamella_path_to_panaroo,
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilli", rownames(Gilliamella_panaroo_out), sep = "_")

Gilliamella_panaroo_core_genes <- rownames(Gilliamella_panaroo_out)[which(Gilliamella_panaroo_out$No..isolates >= 99)]


### First process stringently-competitively mapped breadth data.
Gilliamella_breadth_combined <- read_in_breadth_files(in_path = Gilliamella_path_to_basic_breadth_combined,
                                                      pattern = ".breadth.bedGraph.gz",
                                                      roary_formatted_pangenome = Gilliamella_panaroo_out)

Gilliamella_panaroo_core_genes_consistently_covered <- Gilliamella_panaroo_core_genes[which(Gilliamella_breadth_combined$breadth_summary[Gilliamella_panaroo_core_genes, "min"] >= 0.95)]


# Write out as beds and gffs
Gilliamella_trimmed_panaroo_bed <- read.table("/home/gdouglas/projects/honey_bee/panaroo_pangenomes/clean_panaroo_pangenomes/Gilliamella_panaroo_nonparalogs.100trimmed.bed",
                                              header = FALSE, sep = "\t", stringsAsFactors = FALSE)
rownames(Gilliamella_trimmed_panaroo_bed) <- Gilliamella_trimmed_panaroo_bed$V1

Gilliamella_panaroo_core_genes_consistently_covered_bed <- Gilliamella_trimmed_panaroo_bed[Gilliamella_panaroo_core_genes_consistently_covered, ]
write.table(x = Gilliamella_panaroo_core_genes_consistently_covered_bed,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/consistently_covered/Gilliamella_2019_and_2020_consistent95_core_panaroo.bed",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

Gilliamella_panaroo_core_genes_consistently_covered_gff <- Gilliamella_panaroo_core_genes_consistently_covered_bed
colnames(Gilliamella_panaroo_core_genes_consistently_covered_gff) <- c("scaffold", "start", "stop")
Gilliamella_panaroo_core_genes_consistently_covered_gff$start <- Gilliamella_panaroo_core_genes_consistently_covered_gff$start + 1
Gilliamella_panaroo_core_genes_consistently_covered_gff$creation <- "panaroo"
Gilliamella_panaroo_core_genes_consistently_covered_gff$type <- "CDS"
Gilliamella_panaroo_core_genes_consistently_covered_gff$placeholder1 <- "."
Gilliamella_panaroo_core_genes_consistently_covered_gff$strand <- "+"
Gilliamella_panaroo_core_genes_consistently_covered_gff$placeholder2 <- 0
Gilliamella_panaroo_core_genes_consistently_covered_gff$ID = paste("ID", Gilliamella_panaroo_core_genes_consistently_covered_gff$scaffold, sep = "=")
Gilliamella_panaroo_core_genes_consistently_covered_gff <- Gilliamella_panaroo_core_genes_consistently_covered_gff[, c("scaffold", "creation", "type", "start", "stop",
                                                                                                                                           "placeholder1", "strand", "placeholder2", "ID")]
write.table(x = Gilliamella_panaroo_core_genes_consistently_covered_gff,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/consistently_covered/Gilliamella_2019_and_2020_consistent95_core_panaroo.gff",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")




### Snodgrassella

# Read in reference genome info.
Snodgrassella_panaroo_out <- read.table(Snodgrassella_path_to_panaroo,
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snod", rownames(Snodgrassella_panaroo_out), sep = "_")

Snodgrassella_panaroo_core_genes <- rownames(Snodgrassella_panaroo_out)[which(Snodgrassella_panaroo_out$No..isolates >= 53)]


### First process stringently-competitively mapped breadth data.
Snodgrassella_breadth_combined <- read_in_breadth_files(in_path = Snodgrassella_path_to_basic_breadth_combined,
                                                      pattern = ".breadth.bedGraph.gz",
                                                      roary_formatted_pangenome = Snodgrassella_panaroo_out)

Snodgrassella_panaroo_core_genes_consistently_covered <- Snodgrassella_panaroo_core_genes[which(Snodgrassella_breadth_combined$breadth_summary[Snodgrassella_panaroo_core_genes, "min"] >= 0.75)]


# Write out as beds and gffs
Snodgrassella_trimmed_panaroo_bed <- read.table("/home/gdouglas/projects/honey_bee/panaroo_pangenomes/clean_panaroo_pangenomes/Snodgrassella_panaroo_nonparalogs.100trimmed.bed",
                                              header = FALSE, sep = "\t", stringsAsFactors = FALSE)
rownames(Snodgrassella_trimmed_panaroo_bed) <- Snodgrassella_trimmed_panaroo_bed$V1

Snodgrassella_panaroo_core_genes_consistently_covered_bed <- Snodgrassella_trimmed_panaroo_bed[Snodgrassella_panaroo_core_genes_consistently_covered, ]
write.table(x = Snodgrassella_panaroo_core_genes_consistently_covered_bed,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/consistently_covered/Snodgrassella_2019_and_2020_consistent75_core_panaroo.bed",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

Snodgrassella_panaroo_core_genes_consistently_covered_gff <- Snodgrassella_panaroo_core_genes_consistently_covered_bed
colnames(Snodgrassella_panaroo_core_genes_consistently_covered_gff) <- c("scaffold", "start", "stop")
Snodgrassella_panaroo_core_genes_consistently_covered_gff$start <- Snodgrassella_panaroo_core_genes_consistently_covered_gff$start + 1
Snodgrassella_panaroo_core_genes_consistently_covered_gff$creation <- "panaroo"
Snodgrassella_panaroo_core_genes_consistently_covered_gff$type <- "CDS"
Snodgrassella_panaroo_core_genes_consistently_covered_gff$placeholder1 <- "."
Snodgrassella_panaroo_core_genes_consistently_covered_gff$strand <- "+"
Snodgrassella_panaroo_core_genes_consistently_covered_gff$placeholder2 <- 0
Snodgrassella_panaroo_core_genes_consistently_covered_gff$ID = paste("ID", Snodgrassella_panaroo_core_genes_consistently_covered_gff$scaffold, sep = "=")
Snodgrassella_panaroo_core_genes_consistently_covered_gff <- Snodgrassella_panaroo_core_genes_consistently_covered_gff[, c("scaffold", "creation", "type", "start", "stop",
                                                                                                                       "placeholder1", "strand", "placeholder2", "ID")]
write.table(x = Snodgrassella_panaroo_core_genes_consistently_covered_gff,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/consistently_covered/Snodgrassella_2019_and_2020_consistent75_core_panaroo.gff",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
