rm(list = ls(all.names = TRUE))

library(ape)

# Use mapfile to replace gene tree tip labels so that they contain the genome accessions (which is needed for RANGER-DTL).
raw2clean <- read.table(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/map_raw_to_clean_gene_ids.tsv",
                        sep = "\t", header = FALSE, stringsAsFactors = FALSE, row.names = 1)

replace_tips_labels <- function(in_tree, refound_prefix) {
  
  if (length(grep("refound", in_tree$tip.label)) > 0) {
    in_tree$tip.label[grep("refound", in_tree$tip.label)] <- paste(refound_prefix, in_tree$tip.label[grep("refound", in_tree$tip.label)], sep = "_")
  }
  
  if (length(which(! in_tree$tip.label %in% rownames(raw2clean))) > 0) {
    stop("Tips missing from dataframe!")
  }
  
  in_tree$tip.label <- raw2clean[in_tree$tip.label, "V2"]
  
  return(in_tree)
}

Bifidobacterium_treefiles <- list.files("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/gene_tree_analysis/raw_gene_trees/Bifidobacterium", pattern = ".treefile", full.names = TRUE)
Bifidobacterium_treefiles_basename <- list.files("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/gene_tree_analysis/raw_gene_trees/Bifidobacterium", pattern = ".treefile")
Bifidobacterium_treefiles_basename <- gsub(".treefile", "", Bifidobacterium_treefiles_basename)

for (i in 1:length(Bifidobacterium_treefiles)) {
  
  Bifidobacterium_outfile <- paste("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/gene_tree_analysis/renamed_gene_trees/Bifidobacterium/",
                                    Bifidobacterium_treefiles_basename[i],
                                    ".tree", sep = "")
  
  Bifidobacterium_tree <- replace_tips_labels(in_tree = read.tree(Bifidobacterium_treefiles[i]), refound_prefix = "Bifidobacterium")
   
  write.tree(phy = Bifidobacterium_tree, file = Bifidobacterium_outfile)
  
}



Firm5_treefiles <- list.files("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/gene_tree_analysis/raw_gene_trees/Firm5", pattern = ".treefile", full.names = TRUE)
Firm5_treefiles_basename <- list.files("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/gene_tree_analysis/raw_gene_trees/Firm5", pattern = ".treefile")
Firm5_treefiles_basename <- gsub(".treefile", "", Firm5_treefiles_basename)

for (i in 1:length(Firm5_treefiles)) {
  
  Firm5_outfile <- paste("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/gene_tree_analysis/renamed_gene_trees/Firm5/",
                                   Firm5_treefiles_basename[i],
                                   ".tree", sep = "")
  
  Firm5_tree <- replace_tips_labels(in_tree = read.tree(Firm5_treefiles[i]), refound_prefix = "Firm5")
  
  write.tree(phy = Firm5_tree, file = Firm5_outfile)
  
}




Gilliamella_treefiles <- list.files("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/gene_tree_analysis/raw_gene_trees/Gilliamella", pattern = ".treefile", full.names = TRUE)
Gilliamella_treefiles_basename <- list.files("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/gene_tree_analysis/raw_gene_trees/Gilliamella", pattern = ".treefile")
Gilliamella_treefiles_basename <- gsub(".treefile", "", Gilliamella_treefiles_basename)

for (i in 1:length(Gilliamella_treefiles)) {
  
  Gilliamella_outfile <- paste("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/gene_tree_analysis/renamed_gene_trees/Gilliamella/",
                                   Gilliamella_treefiles_basename[i],
                                   ".tree", sep = "")
  
  Gilliamella_tree <- replace_tips_labels(in_tree = read.tree(Gilliamella_treefiles[i]), refound_prefix = "Gilliamella")
  
  write.tree(phy = Gilliamella_tree, file = Gilliamella_outfile)
  
}




Snodgrassella_treefiles <- list.files("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/gene_tree_analysis/raw_gene_trees/Snodgrassella", pattern = ".treefile", full.names = TRUE)
Snodgrassella_treefiles_basename <- list.files("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/gene_tree_analysis/raw_gene_trees/Snodgrassella", pattern = ".treefile")
Snodgrassella_treefiles_basename <- gsub(".treefile", "", Snodgrassella_treefiles_basename)

for (i in 1:length(Snodgrassella_treefiles)) {
  
  Snodgrassella_outfile <- paste("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/gene_tree_analysis/renamed_gene_trees/Snodgrassella/",
                                   Snodgrassella_treefiles_basename[i],
                                   ".tree", sep = "")
  
  Snodgrassella_tree <- replace_tips_labels(in_tree = read.tree(Snodgrassella_treefiles[i]), refound_prefix = "Snodgrassella")
  
  write.tree(phy = Snodgrassella_tree, file = Snodgrassella_outfile)
  
}
