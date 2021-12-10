rm(list = ls(all.names = TRUE))

library(ape)
library(phytools)

### Perform midpoint rooting on phylotype trees so that they can be used with RANGER-DTL.
### Also simplify the tip names to just be the GCA/GCF ids.

Bifidobacterium_tree <- read.tree("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Bifidobacterium/core_gene_alignment.aln.treefile")
Bifidobacterium_tree_rooted <- midpoint.root(Bifidobacterium_tree)
Bifidobacterium_tree_rooted$tip.label <- gsub("^GCA_", "GCA.", Bifidobacterium_tree_rooted$tip.label)
Bifidobacterium_tree_rooted$tip.label <- gsub("^GCF_", "GCF.", Bifidobacterium_tree_rooted$tip.label)
Bifidobacterium_tree_rooted$tip.label <- gsub("_.*$", "", Bifidobacterium_tree_rooted$tip.label)
write.tree(phy = Bifidobacterium_tree_rooted,
           file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Bifidobacterium/core_gene_alignment.aln.rooted.treefile")


Firm5_tree <- read.tree("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Firm5/core_gene_alignment.aln.treefile")
Firm5_tree_rooted <- midpoint.root(Firm5_tree)
Firm5_tree_rooted$tip.label <- gsub("^GCA_", "GCA.", Firm5_tree_rooted$tip.label)
Firm5_tree_rooted$tip.label <- gsub("^GCF_", "GCF.", Firm5_tree_rooted$tip.label)
Firm5_tree_rooted$tip.label <- gsub("_.*$", "", Firm5_tree_rooted$tip.label)
write.tree(phy = Firm5_tree_rooted, file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Firm5/core_gene_alignment.aln.rooted.treefile")


Gilliamella_tree <- read.tree("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Gilliamella/core_gene_alignment.aln.treefile")
Gilliamella_tree_rooted <- midpoint.root(Gilliamella_tree)
Gilliamella_tree_rooted$tip.label <- gsub("^GCA_", "GCA.", Gilliamella_tree_rooted$tip.label)
Gilliamella_tree_rooted$tip.label <- gsub("^GCF_", "GCF.", Gilliamella_tree_rooted$tip.label)
Gilliamella_tree_rooted$tip.label <- gsub("_.*$", "", Gilliamella_tree_rooted$tip.label)
write.tree(phy = Gilliamella_tree_rooted, file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Gilliamella/core_gene_alignment.aln.rooted.treefile")


Snodgrassella_tree <- read.tree("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Snodgrassella/core_gene_alignment.aln.treefile")
Snodgrassella_tree_rooted <- midpoint.root(Snodgrassella_tree)
Snodgrassella_tree_rooted$tip.label <- gsub("^GCA_", "GCA.", Snodgrassella_tree_rooted$tip.label)
Snodgrassella_tree_rooted$tip.label <- gsub("^GCF_", "GCF.", Snodgrassella_tree_rooted$tip.label)
Snodgrassella_tree_rooted$tip.label <- gsub("_.*$", "", Snodgrassella_tree_rooted$tip.label)
write.tree(phy = Snodgrassella_tree_rooted, file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/Snodgrassella/core_gene_alignment.aln.rooted.treefile")
