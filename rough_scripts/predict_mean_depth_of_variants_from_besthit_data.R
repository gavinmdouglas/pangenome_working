setwd("/home/gdouglas")

pandora_variant_mean_depth <- read.table(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina/pandora_mean_variant_site_depth.tsv",
                                         header = TRUE, sep = "\t", row.names = 1)

phylotypes <- c("Bifidobacterium", "Firm4", "Firm5", "Gilliamella", "Snodgrassella")

for (phylotype in phylotypes) {

  SRR <- "SRR10810031"
  
  besthit_sample_mean_depth <- read.table(paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/map_against_pandora_ref/mean_depth_per_site/",
                                                SRR,
                                                ".mean.bedGraph.gz",
                                                sep = ""),
                                          header = FALSE, sep = "\t", row.names = 1)
  
  combined_df <- data.frame(pandora = pandora_variant_mean_depth[, SRR], besthit = besthit_sample_mean_depth[rownames(pandora_variant_mean_depth), "V4"])
  
  combined_df_subset <- combined_df[-which(combined_df$pandora == 0), ]
  
  variant_mean_depth_lm <- lm(pandora ~ besthit, data = combined_df_subset)

  predict(object = variant_mean_depth_lm, newdata = data.frame(besthit = c(4, 1, 10000)))
  
}