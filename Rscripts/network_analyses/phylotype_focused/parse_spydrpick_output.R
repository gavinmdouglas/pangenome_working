## Parse spydrpick output files and get all significant pairs listed per phylotype in a consistent (and sorted) list saved to an RDS.

rm(list = ls(all.names = TRUE))

setwd("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/")

spydrpick_hits <- list()


Bifidobacterium_spydrpick <- read.table("Bifidobacterium/spydrpick_coevo_results/gene_pa_spydrpick.csv",
                              header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
Bifidobacterium_spydrpick$GeneA <- paste("Bifidobacterium", Bifidobacterium_spydrpick$GeneA, sep = "_")
Bifidobacterium_spydrpick$GeneB <- paste("Bifidobacterium", Bifidobacterium_spydrpick$GeneB, sep = "_")

Bifidobacterium_spydrpick_genes <- Bifidobacterium_spydrpick[, c("GeneA", "GeneB")]
Bifidobacterium_spydrpick_genes_sort <- t(apply(Bifidobacterium_spydrpick_genes, 1, sort))

spydrpick_hits[["Bifidobacterium"]] <- paste(Bifidobacterium_spydrpick_genes_sort[, 1], Bifidobacterium_spydrpick_genes_sort[, 2])


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




# Do the same, but for all non-phylo hits too
rm(list = ls(all.names = TRUE))

setwd("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/")

spydrpick_hits_nonphylo <- list()


Bifidobacterium_spydrpick <- read.table("Bifidobacterium/spydrpick_coevo_results_nophylo/gene_pa_spydrpick.csv",
                              header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
Bifidobacterium_spydrpick$GeneA <- paste("Bifidobacterium", Bifidobacterium_spydrpick$GeneA, sep = "_")
Bifidobacterium_spydrpick$GeneB <- paste("Bifidobacterium", Bifidobacterium_spydrpick$GeneB, sep = "_")

Bifidobacterium_spydrpick_genes <- Bifidobacterium_spydrpick[, c("GeneA", "GeneB")]
Bifidobacterium_spydrpick_genes_sort <- t(apply(Bifidobacterium_spydrpick_genes, 1, sort))

spydrpick_hits_nonphylo[["Bifidobacterium"]] <- paste(Bifidobacterium_spydrpick_genes_sort[, 1], Bifidobacterium_spydrpick_genes_sort[, 2])



Firm4_spydrpick <- read.table("Firm4/spydrpick_coevo_results_nophylo/gene_pa_spydrpick.csv",
                              header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
Firm4_spydrpick$GeneA <- paste("Firm4", Firm4_spydrpick$GeneA, sep = "_")
Firm4_spydrpick$GeneB <- paste("Firm4", Firm4_spydrpick$GeneB, sep = "_")

Firm4_spydrpick_genes <- Firm4_spydrpick[, c("GeneA", "GeneB")]
Firm4_spydrpick_genes_sort <- t(apply(Firm4_spydrpick_genes, 1, sort))

spydrpick_hits_nonphylo[["Firm4"]] <- paste(Firm4_spydrpick_genes_sort[, 1], Firm4_spydrpick_genes_sort[, 2])



Firm5_spydrpick <- read.table("Firm5/spydrpick_coevo_results_nophylo/gene_pa_spydrpick.csv",
                              header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
Firm5_spydrpick$GeneA <- paste("Firm5", Firm5_spydrpick$GeneA, sep = "_")
Firm5_spydrpick$GeneB <- paste("Firm5", Firm5_spydrpick$GeneB, sep = "_")

Firm5_spydrpick_genes <- Firm5_spydrpick[, c("GeneA", "GeneB")]
Firm5_spydrpick_genes_sort <- t(apply(Firm5_spydrpick_genes, 1, sort))

spydrpick_hits_nonphylo[["Firm5"]] <- paste(Firm5_spydrpick_genes_sort[, 1], Firm5_spydrpick_genes_sort[, 2])



Gilliamella_spydrpick <- read.table("Gilliamella/spydrpick_coevo_results_nophylo/gene_pa_spydrpick.csv",
                                    header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
Gilliamella_spydrpick$GeneA <- paste("Gilliamella", Gilliamella_spydrpick$GeneA, sep = "_")
Gilliamella_spydrpick$GeneB <- paste("Gilliamella", Gilliamella_spydrpick$GeneB, sep = "_")

Gilliamella_spydrpick_genes <- Gilliamella_spydrpick[, c("GeneA", "GeneB")]
Gilliamella_spydrpick_genes_sort <- t(apply(Gilliamella_spydrpick_genes, 1, sort))

spydrpick_hits_nonphylo[["Gilliamella"]] <- paste(Gilliamella_spydrpick_genes_sort[, 1], Gilliamella_spydrpick_genes_sort[, 2])



Snodgrassella_spydrpick <- read.table("Snodgrassella/spydrpick_coevo_results_nophylo/gene_pa_spydrpick.csv",
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
Snodgrassella_spydrpick$GeneA <- paste("Snodgrassella", Snodgrassella_spydrpick$GeneA, sep = "_")
Snodgrassella_spydrpick$GeneB <- paste("Snodgrassella", Snodgrassella_spydrpick$GeneB, sep = "_")

Snodgrassella_spydrpick_genes <- Snodgrassella_spydrpick[, c("GeneA", "GeneB")]
Snodgrassella_spydrpick_genes_sort <- t(apply(Snodgrassella_spydrpick_genes, 1, sort))

spydrpick_hits_nonphylo[["Snodgrassella"]] <- paste(Snodgrassella_spydrpick_genes_sort[, 1], Snodgrassella_spydrpick_genes_sort[, 2])


saveRDS(object = spydrpick_hits_nonphylo, file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/spydrpick_hits_nonphylo.rds")
