### Re-ran read mapping with very simple commands to confirm that nothing funky is going on there.

simple_mean_depth <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapping_test/SRR7287194.paired.sort.mean.bedGraph",
                                header = FALSE, sep ="\t", row.names = 1, stringsAsFactors = FALSE)

simple_breadth <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapping_test/SRR7287194.paired.sort.breadth.bedGraph",
                             header = FALSE, sep ="\t", row.names = 1, stringsAsFactors = FALSE)

rownames(simple_mean_depth) <- gsub("Gilliamella", "Gilli", rownames(simple_mean_depth))
rownames(simple_breadth) <- gsub("Gilliamella", "Gilli", rownames(simple_breadth))

rownames(simple_mean_depth) <- gsub("Snodgrassella", "Snod", rownames(simple_mean_depth))
rownames(simple_breadth) <- gsub("Snodgrassella", "Snod", rownames(simple_breadth))


Gilliamella_orig_breadth <- read.table("/data1/gdouglas/projects/honey_bee/Ellegaard_2019_early_work/pangenome_competitive_mapping/bedtools_coverage_out/panaroo_best_hit/Gilliamella_coverage_breadth/SRR7287194.Gilliamella.merged.nonparalog.breadth.bedGraph.gz",
                                       header = FALSE, sep ="\t", row.names = 1, stringsAsFactors = FALSE)

Gilliamella_orig_mean_depth <- read.table("/data1/gdouglas/projects/honey_bee/Ellegaard_2019_early_work/pangenome_competitive_mapping/bedtools_coverage_out/panaroo_best_hit/Gilliamella_mean_depth_per_site/SRR7287194.Gilliamella.merged.nonparalog.mean.bedGraph.gz",
                                          header = FALSE, sep ="\t", row.names = 1, stringsAsFactors = FALSE)

Gilliamella_intersecting <- rownames(Gilliamella_orig_breadth)[which(rownames(Gilliamella_orig_breadth) %in% rownames(simple_mean_depth))]


Snodgrassella_orig_breadth <- read.table("/data1/gdouglas/projects/honey_bee/Ellegaard_2019_early_work/pangenome_competitive_mapping/bedtools_coverage_out/panaroo_best_hit/Snodgrassella_coverage_breadth/SRR7287194.Snodgrassella.merged.nonparalog.breadth.bedGraph.gz",
                                       header = FALSE, sep ="\t", row.names = 1, stringsAsFactors = FALSE)

Snodgrassella_orig_mean_depth <- read.table("/data1/gdouglas/projects/honey_bee/Ellegaard_2019_early_work/pangenome_competitive_mapping/bedtools_coverage_out/panaroo_best_hit/Snodgrassella_mean_depth_per_site/SRR7287194.Snodgrassella.merged.nonparalog.mean.bedGraph.gz",
                                          header = FALSE, sep ="\t", row.names = 1, stringsAsFactors = FALSE)

Snodgrassella_intersecting <- rownames(Snodgrassella_orig_breadth)[which(rownames(Snodgrassella_orig_breadth) %in% rownames(simple_mean_depth))]




plot(simple_mean_depth[Gilliamella_intersecting, "V4"], Gilliamella_orig_mean_depth[Gilliamella_intersecting, "V4"])
plot(simple_breadth[Gilliamella_intersecting, "V7"], Gilliamella_orig_breadth[Gilliamella_intersecting, "V7"])


pandora_output_noncore <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_multisample.filt.noncore.matrix",
                                     header = TRUE, sep = "\t", row.names = 1)

pandora_output <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_multisample.matrix",
                             header = TRUE, sep = "\t", row.names = 1)

rownames(pandora_output) <- gsub(".fa", "", rownames(pandora_output))

pandora_output_core <- pandora_output[which(! rownames(pandora_output) %in% rownames(pandora_output_noncore)), ]

rownames(pandora_output_core) <- gsub("Gilliamella", "Gilli", rownames(pandora_output_core))
rownames(pandora_output_noncore) <- gsub("Gilliamella", "Gilli", rownames(pandora_output_noncore))

rownames(pandora_output_core) <- gsub("Snodgrassella", "Snod", rownames(pandora_output_core))
rownames(pandora_output_noncore) <- gsub("Snodgrassella", "Snod", rownames(pandora_output_noncore))


pandora_and_simple_intersecting_genes_core <- rownames(pandora_output_core)[which(rownames(pandora_output_core) %in% rownames(simple_breadth))]
plot(pandora_output_core[pandora_and_simple_intersecting_genes_core, "SRR7287194"], simple_breadth[pandora_and_simple_intersecting_genes_core, "V7"])

pandora_and_simple_intersecting_genes_noncore <- rownames(pandora_output_noncore)[which(rownames(pandora_output_noncore) %in% rownames(simple_breadth))]
plot(pandora_output_noncore[pandora_and_simple_intersecting_genes_noncore, "SRR7287194"], simple_breadth[pandora_and_simple_intersecting_genes_noncore, "V7"])



pandora_and_Gilliamella_orig_intersecting_genes_core <- rownames(pandora_output_core)[which(rownames(pandora_output_core) %in% rownames(Gilliamella_orig_breadth))]
plot(pandora_output_core[pandora_and_Gilliamella_orig_intersecting_genes_core, "SRR7287194"], Gilliamella_orig_breadth[pandora_and_Gilliamella_orig_intersecting_genes_core, "V7"])


mean(simple_mean_depth[grep("Gilli", pandora_and_simple_intersecting_genes_core, value = TRUE), "V4"])

mean(simple_mean_depth[Gilliamella_intersecting, "V4"])

boxplot(simple_mean_depth[Gilliamella_intersecting, "V4"], Gilliamella_orig_mean_depth[Gilliamella_intersecting, "V4"])

simple_mean_depth[Gilliamella_intersecting, "V4"]


# Do genes with high breadth based on simple approach but that are 0 based on pandora have low depth?

simple_high_breadth <- rownames(simple_breadth[pandora_and_simple_intersecting_genes_core,])[which(simple_breadth[pandora_and_simple_intersecting_genes_core, "V7"] >= 0.9)]

pandora_called <- rownames(pandora_output_core[pandora_and_simple_intersecting_genes_core, ])[which(pandora_output_core[pandora_and_simple_intersecting_genes_core, "SRR7287194"] == 1)]

simple_and_pandora_present <- pandora_called[which(pandora_called %in% simple_high_breadth)]

not_simple_but_pandora_present <- pandora_called[which(! pandora_called %in% simple_high_breadth)]

simple_but_not_pandora_present <- simple_high_breadth[which(! simple_high_breadth %in% pandora_called)]


boxplot(simple_mean_depth[simple_and_pandora_present, "V4"],
        simple_mean_depth[not_simple_but_pandora_present, "V4"],
        simple_mean_depth[simple_but_not_pandora_present, "V4"])





simple_high_breadth <- rownames(simple_breadth[pandora_and_simple_intersecting_genes_noncore,])[which(simple_breadth[pandora_and_simple_intersecting_genes_noncore, "V7"] >= 0.9)]

pandora_called <- rownames(pandora_output_noncore[pandora_and_simple_intersecting_genes_noncore, ])[which(pandora_output_noncore[pandora_and_simple_intersecting_genes_noncore, "SRR7287194"] == 1)]

simple_and_pandora_present <- pandora_called[which(pandora_called %in% simple_high_breadth)]

not_simple_but_pandora_present <- pandora_called[which(! pandora_called %in% simple_high_breadth)]

simple_but_not_pandora_present <- simple_high_breadth[which(! simple_high_breadth %in% pandora_called)]





boxplot(simple_mean_depth[simple_and_pandora_present, "V4"],
        simple_mean_depth[not_simple_but_pandora_present, "V4"],
        simple_mean_depth[simple_but_not_pandora_present, "V4"])




