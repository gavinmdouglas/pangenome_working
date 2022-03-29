rm(list = ls(all.names = TRUE))

indel_breakdown <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/combined_indel_tallies.tsv",
                              header = FALSE, sep = "\t", stringsAsFactors = FALSE)

colnames(indel_breakdown) <- c("gene", "haplotype", "nonframe_insert", "nonframe_del", "frame_insert", "frame_del")

indel_breakdown_by_gene <- aggregate(. ~ gene, data = indel_breakdown, FUN = sum)

panaroo_info <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/combined_panaroo_and_core_genes.rds")

