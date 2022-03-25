rm(list = ls(all.names = TRUE))

# Check to see if dissimilarity-overlap curve analysis on strain presence/abundance results are consistent with universal ecological factors or not.
# Begin by focusing on individual strain datasets at a time, but possibly expand to strains of all species at once, if possible.

# library(devtools)
# install_github("Russel88/DOC")

library("DOC")

strain_abun <- readRDS(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/strain_abun.rds")

Bifidobacterium_asteroides_DOC <- DOC(t(strain_abun$Bifidobacterium_asteroides))
plot(Bifidobacterium_asteroides_DOC$DO$Overlap, Bifidobacterium_asteroides_DOC$DO$rJSD,
     xlab = "Overlap", ylab = "Root Jensen-Shannon Divergence",
     main = "Bifidobacterium_asteroides")

Bombilactobacillus_mellis_DOC <- DOC(t(strain_abun$Bombilactobacillus_mellis))
plot(Bombilactobacillus_mellis_DOC$DO$Overlap, Bombilactobacillus_mellis_DOC$DO$rJSD,
     xlab = "Overlap", ylab = "Root Jensen-Shannon Divergence",
     main = "Bombilactobacillus_mellis")
boxplot(Bombilactobacillus_mellis_DOC$LME$Slope)

Gilliamella_apicola_DOC <- DOC(t(strain_abun$Gilliamella_apicola))
plot(Gilliamella_apicola_DOC$DO$Overlap, Gilliamella_apicola_DOC$DO$rJSD,
     xlab = "Overlap", ylab = "Root Jensen-Shannon Divergence",
     main = "Gilliamella_apicola")
boxplot(Gilliamella_apicola_DOC$LME$Slope)


Lactobacillus_helsingborgensis_DOC <- DOC(t(strain_abun$Lactobacillus_helsingborgensis))
plot(Lactobacillus_helsingborgensis_DOC$DO$Overlap, Lactobacillus_helsingborgensis_DOC$DO$rJSD,
     xlab = "Overlap", ylab = "Root Jensen-Shannon Divergence",
     main = "Lactobacillus_helsingborgensis")
boxplot(Lactobacillus_helsingborgensis_DOC$LME$Slope)


Lactobacillus_melliventris_DOC <- DOC(t(strain_abun$Lactobacillus_melliventris))
plot(Lactobacillus_melliventris_DOC$DO$Overlap, Lactobacillus_melliventris_DOC$DO$rJSD,
     xlab = "Overlap", ylab = "Root Jensen-Shannon Divergence",
     main = "Lactobacillus_melliventris")
boxplot(Lactobacillus_melliventris_DOC$LME$Slope)

