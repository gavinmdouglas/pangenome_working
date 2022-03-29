### Determine majority rule (>=90%) annotation of orthologs based on CDS-level annotations by KEGG.

rm(list = ls(all.names = TRUE))

library(stringr)

setwd("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_ref_only_annot/")

# Read in input
KO_annot <- read.table(file = "KO_microbiota_panaroo_combined_protein_CDS_pro.txt.gz",
                       header = FALSE, sep = "\t", stringsAsFactors = FALSE)

annot2id <- list()
panaroo_out <- list()

for (phylotype in c("Bifidobacterium", "Firm4", "Firm5", "Gilliamella", "Snodgrassella")) {
  
  annot2id[[phylotype]] <- read.table(paste("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/", phylotype, "/gene_data.csv.gz", sep = ""),
                                            header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
  
  rownames(annot2id[[phylotype]]) <- annot2id[[phylotype]]$annotation_id
  
  panaroo_out[[phylotype]] <- read.table(paste("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_out_ref_only/", phylotype, "/gene_presence_absence.csv.gz", sep = ""),
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
  
  rownames(panaroo_out[[phylotype]]) <- paste(phylotype, panaroo_out[[phylotype]]$Gene, sep = "_")
  
}

KO_calls <- list()
for (phylotype in c("Bifidobacterium", "Firm4", "Firm5", "Gilliamella", "Snodgrassella")) {
  
  print(phylotype)
  
  for (ortholog in rownames(panaroo_out[[phylotype]])) {
    ortholog_KO_calls <- c()
    num_CDS <- 0
    for (calls in panaroo_out[[phylotype]][ortholog, 4:ncol(panaroo_out[[phylotype]])]) {
         for (call in str_split(calls, ";")[[1]]) {
          if ((call == "") | is.na(call)) { next }
            num_CDS <- num_CDS + 1
            
            # Remove suffixes indicating short lengths or premature stop codons.
            call <- gsub("_len$", "", call)
            call <- gsub("_stop$", "", call)
            
            # Stop and print name of gene if still not found:
            if (! call %in% rownames(annot2id[[phylotype]])) {
              stop(call) 
            }
            
            call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
            
            if (call_id %in% KO_annot$V1) {
              ortholog_KO_calls <- c(ortholog_KO_calls, KO_annot[which(KO_annot$V1 == call_id), "V2"])
            }
       }
    }
    
    if (length(ortholog_KO_calls) > 0) {
      ortholog_KO_calls_prop <- table(ortholog_KO_calls) / num_CDS
      
      if (length(which(ortholog_KO_calls_prop >= 0.9)) > 0) {
        KO_calls[[ortholog]] <- paste(sort(names(ortholog_KO_calls_prop[which(ortholog_KO_calls_prop >= 0.9)])), collapse = ",")
      }
    }
  }
}

KO_calls_clean <- do.call(c, KO_calls)

KO_calls_out <- data.frame(gene=names(KO_calls_clean), KO=KO_calls_clean)

write.table(x = KO_calls_out,
            file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_ref_only_annot/KO_microbiota_panaroo_orthologs.tsv",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


KO_call_tally_tab <- data.frame(matrix(NA, nrow = 5, ncol = 2))
colnames(KO_call_tally_tab) <- c("count", "percent")
rownames(KO_call_tally_tab) <- c("Bifidobacterium", "Firm4", "Firm5", "Gilliamella", "Snodgrassella")

for (phylotype in rownames(KO_call_tally_tab)) {
  KO_call_tally_tab[phylotype, ] <- c(length(grep(phylotype, KO_calls_out$gene)),
                                        (length(grep(phylotype, KO_calls_out$gene)) / nrow(panaroo_out[[phylotype]])) * 100)
}




##### ##### ##### ##### ##### ##### ##### ##### Sanity checks ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# 
# phylotype <- "Firm4"
# "Firm4_group_1798" 
# panaroo_out[[phylotype]]["Firm4_group_1798", ]
# 
# sample(rownames(panaroo_out$Firm4)[which(! rownames(panaroo_out$Firm4) %in% KO_calls_out$gene)], 3)
# 
# # GEAANGBI_01633              GDAIFKIL_00415              HCAMBJLP_00841
# 
# call <- "GEAANGBI_01633"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# 
# call <- "GDAIFKIL_00415"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# 
# call <- "HCAMBJLP_00841"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# KO_calls_out["Firm4_group_1798",]
# 
# 
# "Firm4_group_2064" 
# panaroo_out[[phylotype]]["Firm4_group_2064", ]
# 
# call <- "GEAANGBI_00840"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# 
# call <- "GDAIFKIL_01569"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# 
# call <- "HCAMBJLP_01634"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]  
# 
# KO_calls_out["Firm4_group_2064",]
# 
# 
# 
# call <- "NEKIOABH_00274"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# KO_calls_out["Firm4_group_1987",]
# 
# 
# 
# panaroo_out[[phylotype]]["Firm4_group_1036", ]
# 
# call <- "GEAANGBI_01530"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# call <- "GDAIFKIL_00514"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# call <- "HCAMBJLP_01057"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# KO_calls_out["Firm4_group_1036", ]
# 
# 
# 
# 
# sample(rownames(panaroo_out$Firm4)[which(! rownames(panaroo_out$Firm4) %in% KO_calls_out$gene)], 3)
# panaroo_out[[phylotype]]["Firm4_group_327", ]
# call <- "PBOALMCJ_00612"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# KO_calls_out["Firm4_group_327", ]
# 
# panaroo_out[[phylotype]]["Firm4_bcgIB", ]
# call <- "HCAMBJLP_01855"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# KO_calls_out["Firm4_bcgIB", ]
# 
# 
# panaroo_out[[phylotype]]["Firm4_yybR", ]
# call <- "GEAANGBI_01028"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# KO_calls_out["Firm4_yybR", ]
# 
# 
# Firm4_frequent_panaroo <- panaroo_out$Firm4[which(rowSums(panaroo_out$Firm4 == "") < 4), ]
# 
# sample(rownames(Firm4_frequent_panaroo)[which(! rownames(Firm4_frequent_panaroo) %in% KO_calls_out$gene)], 3)
# 
# # "Firm4_group_1565"      "Firm4_hepT"            "Firm4_rcsC_2~~~rcsC_1"
# 
# panaroo_out[[phylotype]]["Firm4_group_1565", ]
# # NEKIOABH_00636              GEAANGBI_00384 GDAIFKIL_01085              HCAMBJLP_01264
# call <- "NEKIOABH_00636"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# 
# call <- "GEAANGBI_00384"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# 
# call <- "GDAIFKIL_01085"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# call <- "HCAMBJLP_01264"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# KO_calls_out["Firm4_group_1565", ]
# 
# 
# 
# panaroo_out[[phylotype]]["Firm4_hepT", ]
# # NEKIOABH_01577        GEAANGBI_01509              GDAIFKIL_00539              HCAMBJLP_01079 
# call <- "NEKIOABH_01577"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# 
# call <- "GEAANGBI_01509"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# 
# call <- "GDAIFKIL_00539"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# call <- "HCAMBJLP_01079"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# KO_annot[which(KO_annot$V1 == call_id), "V2"]
# 
# KO_calls_out["Firm4_hepT", ]
# 
# 
# 
# 
# 
# ####### Look into cases where CDS's were annotated as > 1 KO.
# 
# sample(KO_annot$V1[which(duplicated(KO_annot$V1))], 3)
# # "Firm5_16_10_331"     "Firm5_1_0_278"       "Gilliamella_8_54_68"
# 
# phylotype <- "Firm5"
# KO_annot[which(KO_annot$V1 == "Firm5_16_10_331"), ]
# annot_id <- annot2id[[phylotype]][which(annot2id[[phylotype]]$clustering_id == "16_10_331"), "annotation_id"]
# genome_id <- annot2id[[phylotype]][which(annot2id[[phylotype]]$clustering_id == "16_10_331"), "gff_file"]
# 
# CDS_ids <- panaroo_out[[phylotype]][which(panaroo_out[[phylotype]][, genome_id] == annot_id), 4:ncol(panaroo_out[[phylotype]][which(panaroo_out[[phylotype]][, genome_id] == annot_id), ])]
# 
# for (CDS_id in CDS_ids) {
#   call_id <- paste(phylotype, annot2id[[phylotype]][CDS_id, "clustering_id"], sep = "_")
#   print(KO_annot$V2[which(KO_annot$V1 == call_id)])
# }
# 
# KO_calls_out["Firm5_btuD_3~~~btuD_10~~~btuD_8~~~btuD_15~~~btuD_5~~~btuD_11~~~btuD_20", ]
# panaroo_out[[phylotype]]["Firm5_btuD_3~~~btuD_10~~~btuD_8~~~btuD_15~~~btuD_5~~~btuD_11~~~btuD_20", 4:ncol(panaroo_out[[phylotype]])]
# 
# 
# 
# 
# phylotype <- "Firm5"
# KO_annot[which(KO_annot$V1 == "Firm5_1_0_278"), ]
# annot_id <- annot2id[[phylotype]][which(annot2id[[phylotype]]$clustering_id == "1_0_278"), "annotation_id"]
# genome_id <- annot2id[[phylotype]][which(annot2id[[phylotype]]$clustering_id == "1_0_278"), "gff_file"]
# 
# CDS_ids <- panaroo_out[[phylotype]][which(panaroo_out[[phylotype]][, genome_id] == annot_id), 4:ncol(panaroo_out[[phylotype]][which(panaroo_out[[phylotype]][, genome_id] == annot_id), ])]
# 
# gene_id <- paste(phylotype, panaroo_out[[phylotype]][which(panaroo_out[[phylotype]][, genome_id] == annot_id), 1], sep = "_")
# 
# for (CDS_id in CDS_ids) {
#   call_id <- paste(phylotype, annot2id[[phylotype]][CDS_id, "clustering_id"], sep = "_")
#   print(KO_annot$V2[which(KO_annot$V1 == call_id)])
# }
# 
# KO_calls_out[gene_id, ]
# 
# 
# 
# phylotype <- "Gilliamella"
# KO_annot[which(KO_annot$V1 == "Gilliamella_8_54_68"), ]
# annot_id <- annot2id[[phylotype]][which(annot2id[[phylotype]]$clustering_id == "8_54_68"), "annotation_id"]
# genome_id <- annot2id[[phylotype]][which(annot2id[[phylotype]]$clustering_id == "8_54_68"), "gff_file"]
# 
# CDS_ids <- panaroo_out[[phylotype]][which(panaroo_out[[phylotype]][, genome_id] == annot_id), 4:ncol(panaroo_out[[phylotype]][which(panaroo_out[[phylotype]][, genome_id] == annot_id), ])]
# 
# gene_id <- paste(phylotype, panaroo_out[[phylotype]][which(panaroo_out[[phylotype]][, genome_id] == annot_id), 1], sep = "_")
# 
# tmp <- c()
# for (CDS_id in CDS_ids) {
#   call_id <- paste(phylotype, annot2id[[phylotype]][CDS_id, "clustering_id"], sep = "_")
#   #print(KO_annot$V2[which(KO_annot$V1 == call_id)])
#   if(CDS_id == "") { next }
#   tmp <- c(tmp, length(KO_annot$V2[which(KO_annot$V1 == call_id)]))
# }
# 
# KO_calls_out[gene_id, ]
# 
# 
