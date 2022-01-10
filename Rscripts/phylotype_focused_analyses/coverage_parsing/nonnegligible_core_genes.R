### Identify core genes that are consistently mapped to across samples (which could be useful for identifying strains with DESMAN and/or StrainFinder).

rm(list = ls(all.names = TRUE))

setwd("honey_bee_pangenome/data/Ellegaard_panaroo_pangenome_mapping/bowtie2_mapped_gene_coverage/")

source("../../../Rscripts/functions.R")

### Gilliamella

# Read in reference genome info.
Gilliamella_panaroo_out <- read.table("../../binning_output/panaroo_output_dastools_and_ref/Gilliamella_gene_presence_absence_roary-formatted.csv",
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilli", rownames(Gilliamella_panaroo_out), sep = "_")

Gilliamella_panaroo_core_genes <- rownames(Gilliamella_panaroo_out)[which(Gilliamella_panaroo_out$No..isolates >= 99)]


### First process stringently-competitively mapped breadth data.
Gilliamella_breadth_stringent <- read_in_breadth_files(in_path = "stringent_comp_map/Gilliamella_coverage_breadth/",
                                                       pattern = ".breadth.bedGraph.gz",
                                                       roary_formatted_pangenome = Gilliamella_panaroo_out)

Gilliamella_panaroo_core_genes_consistently_covered_stringent <- Gilliamella_panaroo_core_genes[which(Gilliamella_breadth_stringent$breadth_summary[Gilliamella_panaroo_core_genes, "min"] >= 0.95)]


### Then simple approach as well
Gilliamella_breadth_basic <- read_in_breadth_files(in_path = "best_hit_simple_comp_map/Gilliamella_coverage_breadth/",
                                                   pattern = ".breadth.bedGraph.gz",
                                                   roary_formatted_pangenome = Gilliamella_panaroo_out)

Gilliamella_panaroo_core_genes_consistently_covered_basic <- Gilliamella_panaroo_core_genes[which(Gilliamella_breadth_basic$breadth_summary[Gilliamella_panaroo_core_genes, "min"] >= 0.95)]

# Write out as beds and gffs
Gilliamella_trimmed_panaroo_bed <- read.table("../clean_panaroo_pangenomes/Gilliamella_panaroo_nonparalogs.100trimmed.bed.gz",
                                              header = FALSE, sep = "\t", stringsAsFactors = FALSE)
rownames(Gilliamella_trimmed_panaroo_bed) <- Gilliamella_trimmed_panaroo_bed$V1

Gilliamella_panaroo_core_genes_consistently_covered_stringent_bed <- Gilliamella_trimmed_panaroo_bed[Gilliamella_panaroo_core_genes_consistently_covered_stringent, ]
write.table(x = Gilliamella_panaroo_core_genes_consistently_covered_stringent_bed,
            file = "stringent_comp_map/Gilliamella_consistent95_core_panaroo_stringent.bed",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

Gilliamella_panaroo_core_genes_consistently_covered_basic_bed <- Gilliamella_trimmed_panaroo_bed[Gilliamella_panaroo_core_genes_consistently_covered_basic, ]
write.table(x = Gilliamella_panaroo_core_genes_consistently_covered_basic_bed,
            file = "best_hit_simple_comp_map/Gilliamella_consistent95_core_panaroo_basic.bed",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff <- Gilliamella_panaroo_core_genes_consistently_covered_stringent_bed
colnames(Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff) <- c("scaffold", "start", "stop")
Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff$start <- Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff$start + 1
Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff$creation <- "panaroo"
Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff$type <- "CDS"
Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff$placeholder1 <- "."
Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff$strand <- "+"
Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff$placeholder2 <- 0
Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff$ID = paste("ID", Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff$scaffold, sep = "=")
Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff <- Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff[, c("scaffold", "creation", "type", "start", "stop",
                                                                                                                                           "placeholder1", "strand", "placeholder2", "ID")]
write.table(x = Gilliamella_panaroo_core_genes_consistently_covered_stringent_gff,
            file = "stringent_comp_map/Gilliamella_consistent95_core_panaroo_stringent.gff",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

Gilliamella_panaroo_core_genes_consistently_covered_basic_gff <- Gilliamella_panaroo_core_genes_consistently_covered_basic_bed
colnames(Gilliamella_panaroo_core_genes_consistently_covered_basic_gff) <- c("scaffold", "start", "stop")
Gilliamella_panaroo_core_genes_consistently_covered_basic_gff$start <- Gilliamella_panaroo_core_genes_consistently_covered_basic_gff$start + 1
Gilliamella_panaroo_core_genes_consistently_covered_basic_gff$creation <- "panaroo"
Gilliamella_panaroo_core_genes_consistently_covered_basic_gff$type <- "CDS"
Gilliamella_panaroo_core_genes_consistently_covered_basic_gff$placeholder1 <- "."
Gilliamella_panaroo_core_genes_consistently_covered_basic_gff$strand <- "+"
Gilliamella_panaroo_core_genes_consistently_covered_basic_gff$placeholder2 <- 0
Gilliamella_panaroo_core_genes_consistently_covered_basic_gff$ID = paste("ID", Gilliamella_panaroo_core_genes_consistently_covered_basic_gff$scaffold, sep = "=")
Gilliamella_panaroo_core_genes_consistently_covered_basic_gff <- Gilliamella_panaroo_core_genes_consistently_covered_basic_gff[, c("scaffold", "creation", "type", "start", "stop",
                                                                                                                                           "placeholder1", "strand", "placeholder2", "ID")]
write.table(x = Gilliamella_panaroo_core_genes_consistently_covered_basic_gff,
            file = "best_hit_simple_comp_map/Gilliamella_consistent95_core_panaroo_basic.gff",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")


### Snodgrassella

# Read in reference genome info.
Snodgrassella_panaroo_out <- read.table("../../binning_output/panaroo_output_dastools_and_ref/Snodgrassella_gene_presence_absence_roary-formatted.csv",
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snod", rownames(Snodgrassella_panaroo_out), sep = "_")

Snodgrassella_panaroo_core_genes <- rownames(Snodgrassella_panaroo_out)[which(Snodgrassella_panaroo_out$No..isolates >= 53)]


### First process stringently-competitively mapped breadth data.
Snodgrassella_breadth_stringent <- read_in_breadth_files(in_path = "stringent_comp_map/Snodgrassella_coverage_breadth/",
                                                       pattern = ".breadth.bedGraph.gz",
                                                       roary_formatted_pangenome = Snodgrassella_panaroo_out)

Snodgrassella_panaroo_core_genes_consistently_covered_stringent <- Snodgrassella_panaroo_core_genes[which(Snodgrassella_breadth_stringent$breadth_summary[Snodgrassella_panaroo_core_genes, "min"] >= 0.75)]


### Then simple approach as well
Snodgrassella_breadth_basic <- read_in_breadth_files(in_path = "best_hit_simple_comp_map/Snodgrassella_coverage_breadth/",
                                                   pattern = ".breadth.bedGraph.gz",
                                                   roary_formatted_pangenome = Snodgrassella_panaroo_out)

Snodgrassella_panaroo_core_genes_consistently_covered_basic <- Snodgrassella_panaroo_core_genes[which(Snodgrassella_breadth_basic$breadth_summary[Snodgrassella_panaroo_core_genes, "min"] >= 0.75)]

# Write out as bedfiles
Snodgrassella_trimmed_panaroo_bed <- read.table("../clean_panaroo_pangenomes/Snodgrassella_panaroo_nonparalogs.100trimmed.bed.gz",
                                              header = FALSE, sep = "\t", stringsAsFactors = FALSE)
rownames(Snodgrassella_trimmed_panaroo_bed) <- Snodgrassella_trimmed_panaroo_bed$V1

Snodgrassella_panaroo_core_genes_consistently_covered_stringent_bed <- Snodgrassella_trimmed_panaroo_bed[Snodgrassella_panaroo_core_genes_consistently_covered_stringent, ]
write.table(x = Snodgrassella_panaroo_core_genes_consistently_covered_stringent_bed,
            file = "stringent_comp_map/Snodgrassella_consistent75_core_panaroo_stringent.bed",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

Snodgrassella_panaroo_core_genes_consistently_covered_basic_bed <- Snodgrassella_trimmed_panaroo_bed[Snodgrassella_panaroo_core_genes_consistently_covered_basic, ]
write.table(x = Snodgrassella_panaroo_core_genes_consistently_covered_basic_bed,
            file = "best_hit_simple_comp_map/Snodgrassella_consistent75_core_panaroo_basic.bed",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff <- Snodgrassella_panaroo_core_genes_consistently_covered_stringent_bed
colnames(Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff) <- c("scaffold", "start", "stop")
Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff$start <- Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff$start + 1
Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff$creation <- "panaroo"
Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff$type <- "CDS"
Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff$placeholder1 <- "."
Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff$strand <- "+"
Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff$placeholder2 <- 0
Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff$ID = paste("ID", Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff$scaffold, sep = "=")
Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff <- Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff[, c("scaffold", "creation", "type", "start", "stop",
                                                                                                                                           "placeholder1", "strand", "placeholder2", "ID")]
write.table(x = Snodgrassella_panaroo_core_genes_consistently_covered_stringent_gff,
            file = "stringent_comp_map/Snodgrassella_consistent75_core_panaroo_stringent.gff",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff <- Snodgrassella_panaroo_core_genes_consistently_covered_basic_bed
colnames(Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff) <- c("scaffold", "start", "stop")
Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff$start <- Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff$start + 1
Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff$creation <- "panaroo"
Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff$type <- "CDS"
Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff$placeholder1 <- "."
Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff$strand <- "+"
Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff$placeholder2 <- 0
Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff$ID = paste("ID", Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff$scaffold, sep = "=")
Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff <- Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff[, c("scaffold", "creation", "type", "start", "stop",
                                                                                                                                   "placeholder1", "strand", "placeholder2", "ID")]
write.table(x = Snodgrassella_panaroo_core_genes_consistently_covered_basic_gff,
            file = "best_hit_simple_comp_map/Snodgrassella_consistent75_core_panaroo_basic.gff",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
