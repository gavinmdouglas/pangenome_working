rm(list = ls(all.names = TRUE))

library(ggvenn)

spydrpick_hits_phylo <- readRDS(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/spydrpick_hits.rds")

spydrpick_hits_nonphylo <- readRDS(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/spydrpick_hits_nonphylo.rds")

Firm5_hits <- list(Phylo=spydrpick_hits_phylo$Firm5, Nonphylo=spydrpick_hits_nonphylo$Firm5)
ggvenn(Firm5_hits) + ggtitle("Firm5")

Gilliamella_hits <- list(Phylo=spydrpick_hits_phylo$Gilliamella, Nonphylo=spydrpick_hits_nonphylo$Gilliamella)
ggvenn(Gilliamella_hits) + ggtitle("Gilliamella")

Snodgrassella_hits <- list(Phylo=spydrpick_hits_phylo$Snodgrassella, Nonphylo=spydrpick_hits_nonphylo$Snodgrassella)
ggvenn(Snodgrassella_hits) + ggtitle("Snodgrassella")
