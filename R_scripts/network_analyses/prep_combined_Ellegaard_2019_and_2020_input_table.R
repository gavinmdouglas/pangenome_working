rm(list = ls(all.names = TRUE))

pandora_output <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present//pandora_multisample.matrix",
                             header = TRUE, sep = "\t", row.names = 1)

panaroo_only_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/panaroo_only_potential_core.rds")

checkm_panaroo_passed_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/checkm_panaroo_passed_potential_core.rds")

core_genes <- c(panaroo_only_core, checkm_panaroo_passed_core)
all_core_genes <- do.call(c, core_genes)
names(all_core_genes) <- NULL

pandora_output_noncore <- pandora_output[-which(rownames(pandora_output) %in% all_core_genes), ]

write.table(x = pandora_output_noncore,
            file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina_all_present/pandora_multisample.noncore.matrix",
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

