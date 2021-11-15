### Remove samples with insufficient coverage per genus from the input file for pandora.

rm(list = ls(all.names = TRUE))

pandora_output_filt_illumina <- read.table(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina/pandora_multisample.filt.matrix",
                                            sep = "\t", row.names = 1, header = TRUE)

pandora_output_filt_default <- read.table(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_orig/pandora_multisample.filt.matrix",
                                           sep = "\t", row.names = 1, header = TRUE)

intersecting_genes <- rownames(pandora_output_filt_illumina)[which(rownames(pandora_output_filt_illumina) %in% rownames(pandora_output_filt_default))]

num_same <- c()

for (s in colnames(pandora_output_filt_illumina)) {
  num_same <- c(num_same, length(which(pandora_output_filt_illumina[intersecting_genes, s] == pandora_output_filt_default[intersecting_genes, s])))
}

min(num_same / length(intersecting_genes))

# Mean of 0.9982 and min of 0.9956. Seems to have not made much of a difference!
