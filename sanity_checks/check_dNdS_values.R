rm(list = ls(all.names = TRUE))

# Check a few dN/dS values with ape compared with my script.

library(ape)

test1 <- read.FASTA("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/sf_prepped_genes_codon_alignments/Snodgrassella_alvi/Snodgrassella_alvi_yqjC.fa", type = "DNA")

dnds(test1, quiet = TRUE)

