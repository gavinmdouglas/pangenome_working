rm(list = ls(all.names = TRUE))

source("honey_bee_pangenome/Rscripts/functions.R")

# Read in reference genome info.
Gilliamella_panaroo_out <- read.table("honey_bee_pangenome/data/binning_output/panaroo_output_dastools_and_ref/Gilliamella_gene_presence_absence_roary-formatted.csv",
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Gilliamella_panaroo_out) <- paste("Gilli", rownames(Gilliamella_panaroo_out), sep = "_")

Gilliamella_basic_breadth <- read_in_breadth_files(in_path = "honey_bee_pangenome/data/Ellegaard_panaroo_pangenome_mapping/bowtie2_mapped_gene_coverage/best_hit_simple_comp_map/Gilliamella_coverage_breadth/",
                                                   pattern = ".breadth.bedGraph.gz",
                                                   roary_formatted_pangenome = Gilliamella_panaroo_out)

Gilliamella_no.comp.map_breadth <- read_in_breadth_files(in_path = "honey_bee_pangenome/data/tests/test_coverage_breadth_no_competitive_mapping/Gilliamella_only_map_out_coverage_breadth//",
                                                           pattern = ".breadth.bedGraph.gz",
                                                           roary_formatted_pangenome = Gilliamella_panaroo_out)

Gilliamella_basic_panaroo_core <- rownames(Gilliamella_panaroo_out)[which(Gilliamella_panaroo_out$No..isolates >= 99)]

rownames(Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287214) <- Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287214$gene
rownames(Gilliamella_basic_breadth$breadth_by_sample$SRR7287214) <- Gilliamella_basic_breadth$breadth_by_sample$SRR7287214$gene

rownames(Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287220) <- Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287220$gene
rownames(Gilliamella_basic_breadth$breadth_by_sample$SRR7287220) <- Gilliamella_basic_breadth$breadth_by_sample$SRR7287220$gene

rownames(Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287234) <- Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287234$gene
rownames(Gilliamella_basic_breadth$breadth_by_sample$SRR7287234) <- Gilliamella_basic_breadth$breadth_by_sample$SRR7287234$gene


rownames(Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287214) <- Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287214$gene
rownames(Gilliamella_basic_breadth$breadth_by_sample$SRR7287214) <- Gilliamella_basic_breadth$breadth_by_sample$SRR7287214$gene


boxplot(Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287214[Gilliamella_basic_panaroo_core, "V7"],
        Gilliamella_basic_breadth$breadth_by_sample$SRR7287214[Gilliamella_basic_panaroo_core, "V7"])


plot(Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287214[Gilliamella_basic_panaroo_core, "V7"],
        Gilliamella_basic_breadth$breadth_by_sample$SRR7287214[Gilliamella_basic_panaroo_core, "V7"])


which(Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287214[Gilliamella_basic_panaroo_core, "V7"] -
        Gilliamella_basic_breadth$breadth_by_sample$SRR7287214[Gilliamella_basic_panaroo_core, "V7"] == min(Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287214[Gilliamella_basic_panaroo_core, "V7"]-
                                                                                                            Gilliamella_basic_breadth$breadth_by_sample$SRR7287214[Gilliamella_basic_panaroo_core, "V7"], na.rm = TRUE))



plot(Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287234[Gilliamella_basic_panaroo_core, "V7"],
     Gilliamella_basic_breadth$breadth_by_sample$SRR7287234[Gilliamella_basic_panaroo_core, "V7"])

plot(Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287234[Gilliamella_basic_panaroo_core, "V7"],
     Gilliamella_basic_breadth$breadth_by_sample$SRR7287234[Gilliamella_basic_panaroo_core, "V7"])








# Read in reference genome info.
Snodgrassella_panaroo_out <- read.table("honey_bee_pangenome/data/binning_output/panaroo_output_dastools_and_ref/Snodgrassella_gene_presence_absence_roary-formatted.csv",
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)
rownames(Snodgrassella_panaroo_out) <- paste("Snod", rownames(Snodgrassella_panaroo_out), sep = "_")

Snodgrassella_basic_breadth <- read_in_breadth_files(in_path = "honey_bee_pangenome/data/Ellegaard_panaroo_pangenome_mapping/bowtie2_mapped_gene_coverage/best_hit_simple_comp_map/Snodgrassella_coverage_breadth/",
                                                   pattern = ".breadth.bedGraph.gz",
                                                   roary_formatted_pangenome = Snodgrassella_panaroo_out)

Snodgrassella_no.comp.map_breadth <- read_in_breadth_files(in_path = "honey_bee_pangenome/data/tests/test_coverage_breadth_no_competitive_mapping/Snodgrassella_only_map_out_coverage_breadth//",
                                                         pattern = ".breadth.bedGraph.gz",
                                                         roary_formatted_pangenome = Snodgrassella_panaroo_out)

Snodgrassella_basic_panaroo_core <- rownames(Snodgrassella_panaroo_out)[which(Snodgrassella_panaroo_out$No..isolates >= 53)]

rownames(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287214) <- Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287214$gene
rownames(Snodgrassella_basic_breadth$breadth_by_sample$SRR7287214) <- Snodgrassella_basic_breadth$breadth_by_sample$SRR7287214$gene

rownames(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287220) <- Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287220$gene
rownames(Snodgrassella_basic_breadth$breadth_by_sample$SRR7287220) <- Snodgrassella_basic_breadth$breadth_by_sample$SRR7287220$gene

rownames(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287234) <- Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287234$gene
rownames(Snodgrassella_basic_breadth$breadth_by_sample$SRR7287234) <- Snodgrassella_basic_breadth$breadth_by_sample$SRR7287234$gene


rownames(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287214) <- Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287214$gene
rownames(Snodgrassella_basic_breadth$breadth_by_sample$SRR7287214) <- Snodgrassella_basic_breadth$breadth_by_sample$SRR7287214$gene


boxplot(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287214[Snodgrassella_basic_panaroo_core, "V7"],
        Snodgrassella_basic_breadth$breadth_by_sample$SRR7287214[Snodgrassella_basic_panaroo_core, "V7"])


plot(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287214[Snodgrassella_basic_panaroo_core, "V7"],
     Snodgrassella_basic_breadth$breadth_by_sample$SRR7287214[Snodgrassella_basic_panaroo_core, "V7"])


which(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287214[Snodgrassella_basic_panaroo_core, "V7"] -
        Snodgrassella_basic_breadth$breadth_by_sample$SRR7287214[Snodgrassella_basic_panaroo_core, "V7"] == min(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287214[Snodgrassella_basic_panaroo_core, "V7"]-
                                                                                                              Snodgrassella_basic_breadth$breadth_by_sample$SRR7287214[Snodgrassella_basic_panaroo_core, "V7"], na.rm = TRUE))



plot(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287234[Snodgrassella_basic_panaroo_core, "V7"],
     Snodgrassella_basic_breadth$breadth_by_sample$SRR7287234[Snodgrassella_basic_panaroo_core, "V7"])

plot(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287234[Snodgrassella_basic_panaroo_core, "V7"],
     Snodgrassella_basic_breadth$breadth_by_sample$SRR7287234[Snodgrassella_basic_panaroo_core, "V7"])







### Make clean plot

par(mfrow = c(2, 3))

plot(Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287214[Gilliamella_basic_panaroo_core, "V7"],
     Gilliamella_basic_breadth$breadth_by_sample$SRR7287214[Gilliamella_basic_panaroo_core, "V7"],
     main = "Gilliamella - SRR7287214",
     xlab = "No comp. map",
     ylab = "Comp. map")

plot(Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287220[Gilliamella_basic_panaroo_core, "V7"],
     Gilliamella_basic_breadth$breadth_by_sample$SRR7287220[Gilliamella_basic_panaroo_core, "V7"],
     main = "Gilliamella - SRR7287220",
     xlab = "No comp. map",
     ylab = "Comp. map")

plot(Gilliamella_no.comp.map_breadth$breadth_by_sample$SRR7287234[Gilliamella_basic_panaroo_core, "V7"],
     Gilliamella_basic_breadth$breadth_by_sample$SRR7287234[Gilliamella_basic_panaroo_core, "V7"],
     main = "Gilliamella - SRR7287234",
     xlab = "No comp. map",
     ylab = "Comp. map")


plot(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287214[Snodgrassella_basic_panaroo_core, "V7"],
     Snodgrassella_basic_breadth$breadth_by_sample$SRR7287214[Snodgrassella_basic_panaroo_core, "V7"],
     main = "Snodgrassella - SRR7287214",
     xlab = "No comp. map",
     ylab = "Comp. map")

plot(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287220[Snodgrassella_basic_panaroo_core, "V7"],
     Snodgrassella_basic_breadth$breadth_by_sample$SRR7287220[Snodgrassella_basic_panaroo_core, "V7"],
     main = "Snodgrassella - SRR7287220",
     xlab = "No comp. map",
     ylab = "Comp. map")

plot(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287234[Snodgrassella_basic_panaroo_core, "V7"],
     Snodgrassella_basic_breadth$breadth_by_sample$SRR7287234[Snodgrassella_basic_panaroo_core, "V7"],
     main = "Snodgrassella - SRR7287234",
     xlab = "No comp. map",
     ylab = "Comp. map")

par(mfrow = c(1, 1))




plot(Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287220[Snodgrassella_basic_panaroo_core, "V7"],
     Snodgrassella_no.comp.map_breadth$breadth_by_sample$SRR7287234[Snodgrassella_basic_panaroo_core, "V7"])