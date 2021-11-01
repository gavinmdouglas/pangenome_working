## Parse spydrpick output files and get all significant pairs listed per phylotype in a consistent (and sorted) list saved to an RDS.

rm(list = ls(all.names = TRUE))

setwd("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/")

spydrpick_hits <- list()


Firm5_spydrpick <- read.table("Firm5/spydrpick_coevo_results/gene_pa_spydrpick.csv",
                              header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
Firm5_spydrpick$GeneA <- paste("Firm5", Firm5_spydrpick$GeneA, sep = "_")
Firm5_spydrpick$GeneB <- paste("Firm5", Firm5_spydrpick$GeneB, sep = "_")

Firm5_spydrpick_genes <- Firm5_spydrpick[, c("GeneA", "GeneB")]
Firm5_spydrpick_genes_sort <- t(apply(Firm5_spydrpick_genes, 1, sort))

spydrpick_hits[["Firm5"]] <- paste(Firm5_spydrpick_genes_sort[, 1], Firm5_spydrpick_genes_sort[, 2])



Gilliamella_spydrpick <- read.table("Gilliamella/spydrpick_coevo_results/gene_pa_spydrpick.csv",
                              header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
Gilliamella_spydrpick$GeneA <- paste("Gilliamella", Gilliamella_spydrpick$GeneA, sep = "_")
Gilliamella_spydrpick$GeneB <- paste("Gilliamella", Gilliamella_spydrpick$GeneB, sep = "_")

Gilliamella_spydrpick_genes <- Gilliamella_spydrpick[, c("GeneA", "GeneB")]
Gilliamella_spydrpick_genes_sort <- t(apply(Gilliamella_spydrpick_genes, 1, sort))

spydrpick_hits[["Gilliamella"]] <- paste(Gilliamella_spydrpick_genes_sort[, 1], Gilliamella_spydrpick_genes_sort[, 2])



Snodgrassella_spydrpick <- read.table("Snodgrassella/spydrpick_coevo_results/gene_pa_spydrpick.csv",
                              header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
Snodgrassella_spydrpick$GeneA <- paste("Snodgrassella", Snodgrassella_spydrpick$GeneA, sep = "_")
Snodgrassella_spydrpick$GeneB <- paste("Snodgrassella", Snodgrassella_spydrpick$GeneB, sep = "_")

Snodgrassella_spydrpick_genes <- Snodgrassella_spydrpick[, c("GeneA", "GeneB")]
Snodgrassella_spydrpick_genes_sort <- t(apply(Snodgrassella_spydrpick_genes, 1, sort))

spydrpick_hits[["Snodgrassella"]] <- paste(Snodgrassella_spydrpick_genes_sort[, 1], Snodgrassella_spydrpick_genes_sort[, 2])



saveRDS(object = spydrpick_hits, file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/spydrpick_hits.rds")
