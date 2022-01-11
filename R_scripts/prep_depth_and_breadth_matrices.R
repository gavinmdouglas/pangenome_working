### Prep presence absence and depth matrices for both Ellegaard 2020 and 2019 samples

rm(list = ls(all.names = TRUE))

library("plyr")

source("/home/gdouglas/scripts/pangenome_working/Rscripts/functions.R")

# First get processed breadth of coverage table.
breadth_outdir <- "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/map_against_pandora_ref/coverage_breadth/"

combined_panaroo_and_core_genes <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/combined_panaroo_and_core_genes.rds")

breadth_summaries <- read_in_breadth_files(in_path = breadth_outdir,
                                           pattern = ".breadth.bedGraph.gz",
                                           roary_formatted_pangenome = combined_panaroo_and_core_genes$panaroo_all_phylotypes)

breadth_by_sample <- list()

all_samples <- names(breadth_summaries$breadth_by_sample)

for (s in all_samples) {
  breadth_by_sample[[s]] <- breadth_summaries$breadth_by_sample[[s]][, c("gene", "V7"), drop = FALSE]
  colnames(breadth_by_sample[[s]]) <- c("gene", s)
}

all_breadth <- join_all(breadth_by_sample, by="gene", type='left')
rownames(all_breadth) <- all_breadth$gene
all_breadth <- all_breadth[, -which(colnames(all_breadth) == "gene")]

saveRDS(object = breadth_summaries,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/map_against_pandora_ref/clean_tables/Ellegaard.2019.2020_pandora_panaroo_coverage_breadth_summaries.rds")

saveRDS(object = all_breadth,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/map_against_pandora_ref/clean_tables/Ellegaard.2019.2020_pandora_panaroo_coverage_breadth.rds")



# Then get processed mean depth table.
depth_outdir <- "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/map_against_pandora_ref/mean_depth_per_site/"

depth_summaries <- read_in_depth_files(in_path = depth_outdir,
                                       pattern = ".mean.bedGraph.gz",
                                       roary_formatted_pangenome = combined_panaroo_and_core_genes$panaroo_all_phylotypes)

depth_by_sample <- list()

for (s in all_samples) {
  depth_by_sample[[s]] <- depth_summaries$depth_by_sample[[s]][, c("gene", "V4"), drop = FALSE]
  colnames(depth_by_sample[[s]]) <- c("gene", s)
}


all_depth <- join_all(depth_by_sample, by="gene", type='left')
rownames(all_depth) <- all_depth$gene
all_depth <- all_depth[, -which(colnames(all_depth) == "gene")]

saveRDS(object = depth_summaries,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/map_against_pandora_ref/clean_tables/Ellegaard.2019.2020_pandora_panaroo_mean_depth_summaries.rds")

saveRDS(object = all_depth,
        file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/map_against_pandora_ref/clean_tables/Ellegaard.2019.2020_pandora_panaroo_mean_depth.rds")
