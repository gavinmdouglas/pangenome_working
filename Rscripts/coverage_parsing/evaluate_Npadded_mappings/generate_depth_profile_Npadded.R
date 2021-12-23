### For one test file, check to see to what degree and in what position N bases have been mapped to.
### This comparison will be relative to all other bases in the locus.

rm(list = ls(all.names = TRUE))

setwd("honey_bee_pangenome/data/Ellegaard_Carrie_pangenome_mapping/Nfilled_mapping_test/")

Gilliamella_pangenome_annot <- read.table("Gilli_annot_Nfilled_bamreadcount_format.tsv.gz",
                                          sep = "\t", stringsAsFactors = FALSE)
colnames(Gilliamella_pangenome_annot) <- c("subfam", "start", "stop")

Gilliamella_bam.readcount <- read.table("SRR7287226.merged.bam-readcount.simple.Nfilled.tsv.gz",
                                        sep = "\t", stringsAsFactors = FALSE)
colnames(Gilliamella_bam.readcount) <- c("subfam", "position", "ref_base", "depth")

# Only keep subfams in annot that are in bam.readcount
Gilliamella_pangenome_annot <- Gilliamella_pangenome_annot[which(Gilliamella_pangenome_annot$subfam %in% Gilliamella_bam.readcount$subfam), ]

# Look at distribution of max depth per subfamily.
Gilliamella_subfam_max_depth <- c()

for (Gilliamella_subfam_id in Gilliamella_pangenome_annot$subfam) {
  Gilliamella_subfam_max_depth <- c(Gilliamella_subfam_max_depth, max(Gilliamella_bam.readcount[which(Gilliamella_bam.readcount$subfam == Gilliamella_subfam_id), "depth"]))
}

Gilliamella_subfams2keep <- Gilliamella_pangenome_annot$subfam[which(Gilliamella_subfam_max_depth >= 20)]

# Need to make sure that all loci used for this analysis are >= 500 bp
rownames(Gilliamella_pangenome_annot) <- Gilliamella_pangenome_annot$subfam
Gilliamella_pangenome_annot <- Gilliamella_pangenome_annot[Gilliamella_subfams2keep, ]
Gilliamella_pangenome_annot <- Gilliamella_pangenome_annot[which(Gilliamella_pangenome_annot$stop >= 500), ]
Gilliamella_subfams2keep <- Gilliamella_pangenome_annot$subfam
  
# Save object of tested subfamilies in order to keep this consistent with 
saveRDS(object = Gilliamella_subfams2keep, file = "Gilliamella_subfams_tested.rds")

leading_N_labels <- 1:50
leading_actual_labels <- 51:300
trailing_actual_labels <- -300:-51
trailing_N_labels <- -50:-1

xlabels <- c(leading_N_labels, leading_actual_labels, trailing_actual_labels, trailing_N_labels)
Gilliamella_depth_by_pos <- data.frame(matrix(NA, nrow = length(xlabels), ncol = 3))

colnames(Gilliamella_depth_by_pos) <- c("xlabels", "raw_depth", "num_instances")
rownames(Gilliamella_depth_by_pos) <- as.character(xlabels)
Gilliamella_depth_by_pos$xlabels <- xlabels
Gilliamella_depth_by_pos$raw_depth <- 0
Gilliamella_depth_by_pos$num_instances <- 0

for (Gilliamella_subfam_id in Gilliamella_subfams2keep) {

  Gilliamella_locus_length <- as.numeric(Gilliamella_pangenome_annot[which(Gilliamella_pangenome_annot$subfam == Gilliamella_subfam_id), "stop"])
  
  Gilliamella_leading_N_i <- 1:50
  Gilliamella_leading_actual_i <- 51:300

  Gilliamella_trailing_actual_i <- (Gilliamella_locus_length - 299):(Gilliamella_locus_length - 50)
  Gilliamella_trailing_N_i <- (Gilliamella_locus_length - 49):Gilliamella_locus_length
  
  Gilliamella_bam.readcount_subset <- Gilliamella_bam.readcount[which(Gilliamella_bam.readcount$subfam == Gilliamella_subfam_id), ]
  
  Gilliamella_leading_N_i_intersect <- Gilliamella_leading_N_i[which(Gilliamella_leading_N_i %in% Gilliamella_bam.readcount_subset$position)]
  if (length(Gilliamella_leading_N_i_intersect) > 0) {
    for (leading_N_i in Gilliamella_leading_N_i_intersect) {
      Gilliamella_bam.readcount_subset_leading_N_i <- Gilliamella_bam.readcount_subset[which(Gilliamella_bam.readcount_subset$position == leading_N_i), ]
      Gilliamella_tmp_leading_N_label <- as.character(leading_N_labels[which(Gilliamella_leading_N_i %in% leading_N_i)])
      Gilliamella_depth_by_pos[Gilliamella_tmp_leading_N_label, c("raw_depth", "num_instances")] <- Gilliamella_depth_by_pos[Gilliamella_tmp_leading_N_label, c("raw_depth", "num_instances")] + c(Gilliamella_bam.readcount_subset_leading_N_i$depth, 1)
    }
  }
 
  Gilliamella_leading_actual_i_intersect <- Gilliamella_leading_actual_i[which(Gilliamella_leading_actual_i %in% Gilliamella_bam.readcount_subset$position)]
  if (length(Gilliamella_leading_actual_i_intersect) > 0) {
    for (leading_actual_i in Gilliamella_leading_actual_i_intersect) {
      Gilliamella_bam.readcount_subset_leading_actual_i <- Gilliamella_bam.readcount_subset[which(Gilliamella_bam.readcount_subset$position == leading_actual_i), ]
      Gilliamella_tmp_leading_actual_label <- as.character(leading_actual_labels[which(Gilliamella_leading_actual_i %in% leading_actual_i)])
      Gilliamella_depth_by_pos[Gilliamella_tmp_leading_actual_label, c("raw_depth", "num_instances")] <- Gilliamella_depth_by_pos[Gilliamella_tmp_leading_actual_label, c("raw_depth", "num_instances")] + c(Gilliamella_bam.readcount_subset_leading_actual_i$depth, 1)
    }
  }
  
  Gilliamella_trailing_actual_i_intersect <- Gilliamella_trailing_actual_i[which(Gilliamella_trailing_actual_i %in% Gilliamella_bam.readcount_subset$position)]
  if (length(Gilliamella_trailing_actual_i_intersect) > 0) {
    for (trailing_actual_i in Gilliamella_trailing_actual_i_intersect) {
      Gilliamella_bam.readcount_subset_trailing_actual_i <- Gilliamella_bam.readcount_subset[which(Gilliamella_bam.readcount_subset$position == trailing_actual_i), ]
      Gilliamella_tmp_trailing_actual_label <- as.character(trailing_actual_labels[which(Gilliamella_trailing_actual_i %in% trailing_actual_i)])
      Gilliamella_depth_by_pos[Gilliamella_tmp_trailing_actual_label, c("raw_depth", "num_instances")] <- Gilliamella_depth_by_pos[Gilliamella_tmp_trailing_actual_label, c("raw_depth", "num_instances")] + c(Gilliamella_bam.readcount_subset_trailing_actual_i$depth, 1)
    }
  }
  
  Gilliamella_trailing_N_i_intersect <- Gilliamella_trailing_N_i[which(Gilliamella_trailing_N_i %in% Gilliamella_bam.readcount_subset$position)]
  if (length(Gilliamella_trailing_N_i_intersect) > 0) {
    for (trailing_N_i in Gilliamella_trailing_N_i_intersect) {
      Gilliamella_bam.readcount_subset_trailing_N_i <- Gilliamella_bam.readcount_subset[which(Gilliamella_bam.readcount_subset$position == trailing_N_i), ]
      Gilliamella_tmp_trailing_N_label <- as.character(trailing_N_labels[which(Gilliamella_trailing_N_i %in% trailing_N_i)])
      Gilliamella_depth_by_pos[Gilliamella_tmp_trailing_N_label, c("raw_depth", "num_instances")] <- Gilliamella_depth_by_pos[Gilliamella_tmp_trailing_N_label, c("raw_depth", "num_instances")] + c(Gilliamella_bam.readcount_subset_trailing_N_i$depth, 1)
    }
  }
}

saveRDS(object = Gilliamella_depth_by_pos,
        file = "Nfilled_Gilliamella_depth_by_pos.rds")


