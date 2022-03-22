# Read in all SF output and save as RDS, because it takes a long time to do this step.

rm(list = ls(all.names = TRUE))

species_present <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species_present.txt", header = FALSE, stringsAsFactors = FALSE)$V1

prepped_genes <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/all_genes_prepped.txt",
                            header = FALSE, stringsAsFactors = FALSE)$V1

invariant_genes <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/invariant_genes.txt",
                              header = FALSE, stringsAsFactors = FALSE)$V1

presence_matrix <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.matrix",
                              header = TRUE, sep = "\t", row.names = 1)


outfolder <- "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_output"

AIC <- list()
OTU <- list()

for (species in species_present) {
  AIC[[species]] <- list()
  OTU[[species]] <- list()
}

for (gene in prepped_genes) {
  
  gene_split <- strsplit(x = gene, split = "_")[[1]]
  
  species <- paste(gene_split[1], gene_split[2], sep = "_")
  
  if (species == "Bifidobacterium_coryneforme") {
    species <- "Bifidobacterium_coryneforme_indicum" 
  }
  
  gene_sf_output <- list.files(paste(outfolder, species, gene, sep = "/"), full.names = TRUE)
  
  AIC_file <- grep("strain_fit_summary.tsv", gene_sf_output, value = TRUE)
  
  OTU_file <- grep("otu_table", gene_sf_output, value = TRUE)
  
  sample_file <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/prepped_input/", gene, "_samples.txt", sep = "")
  
  AIC[[species]][[gene]] <- read.table(file = AIC_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  OTU[[species]][[gene]] <- read.table(file = OTU_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  rownames(OTU[[species]][[gene]]) <- read.table(file = sample_file, header = FALSE, stringsAsFactors = FALSE)$V1
  
  OTU[[species]][[gene]][OTU[[species]][[gene]] < 0.01] <- 0
  
  OTU[[species]][[gene]] <- data.frame(sweep(OTU[[species]][[gene]], 1, rowSums(OTU[[species]][[gene]]), '/')) * 100

}

saveRDS(object = AIC, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/all_genes_AIC.rds")
saveRDS(object = OTU, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/output_summaries/all_genes_OTU.rds")

