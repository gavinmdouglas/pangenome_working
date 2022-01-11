### Now that consensus CAZy annotations have been parsed for each individual CDS sequences,
### now can focus on getting annotations for overall orthologs for which groups of CDS sequences belong to.

rm(list = ls(all.names = TRUE))

library(stringr)

setwd("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_ref_only_annot/")

# Read in input
dbCAN2_all_hits_atleast2 <- read.table(file = "dbCAN2_microbiota_panaroo_combined_protein_CDS_consensus.tsv.gz",
                                       header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

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

CAZy_calls <- list()
for (phylotype in c("Bifidobacterium", "Firm4", "Firm5", "Gilliamella", "Snodgrassella")) {
  
  print(phylotype)
  
  for (ortholog in rownames(panaroo_out[[phylotype]])) {
    ortholog_CAZy_calls <- c()
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
            
            if (call_id %in% rownames(dbCAN2_all_hits_atleast2)) {
              ortholog_CAZy_calls <- c(ortholog_CAZy_calls, str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]])
            }
       }
    }
    
    if (length(ortholog_CAZy_calls) > 0) {
      ortholog_CAZy_calls_prop <- table(ortholog_CAZy_calls) / num_CDS
      
      if (length(which(ortholog_CAZy_calls_prop >= 0.9)) > 0) {
        CAZy_calls[[ortholog]] <- paste(sort(names(ortholog_CAZy_calls_prop[which(ortholog_CAZy_calls_prop >= 0.9)])), collapse = ",")
      }
    }
  }
}

CAZy_calls_clean <- do.call(c, CAZy_calls)

CAZy_calls_out <- data.frame(gene=names(CAZy_calls_clean), CAZy=CAZy_calls_clean)

write.table(x = CAZy_calls_out,
            file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/panaroo_ref_only_annot/CAZy_microbiota_panaroo_orthologs.tsv",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


CAZy_call_tally_tab <- data.frame(matrix(NA, nrow = 5, ncol = 2))
colnames(CAZy_call_tally_tab) <- c("count", "percent")
rownames(CAZy_call_tally_tab) <- c("Bifidobacterium", "Firm4", "Firm5", "Gilliamella", "Snodgrassella")

for (phylotype in rownames(CAZy_call_tally_tab)) {
  CAZy_call_tally_tab[phylotype, ] <- c(length(grep(phylotype, CAZy_calls_out$gene)),
                                        (length(grep(phylotype, CAZy_calls_out$gene)) / nrow(panaroo_out[[phylotype]])) * 100)
}




##### ##### ##### ##### ##### ##### ##### ##### Sanity checks ##### ##### ##### ##### ##### ##### ##### ##### ##### 

# phylotype <- "Firm4"
# 
# "Firm4_group_1798" 
# panaroo_out[[phylotype]]["Firm4_group_1798", ]
# 
# # GEAANGBI_01633              GDAIFKIL_00415              HCAMBJLP_00841
# 
# call <- "GEAANGBI_01633"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# 
# call <- "GDAIFKIL_00415"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# 
# call <- "HCAMBJLP_00841"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# CAZy_calls["Firm4_group_1798"]
# 
# 
# 
# "Firm4_group_2064" 
# panaroo_out[[phylotype]]["Firm4_group_2064", ]
# 
# call <- "GEAANGBI_00840"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# 
# call <- "GDAIFKIL_01569"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# 
# call <- "HCAMBJLP_01634"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]  
# CAZy_calls["Firm4_group_2064"]
# 
# 
# 
# panaroo_out[[phylotype]]["Firm4_group_1987", ]
# 
# call <- "NEKIOABH_00274"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# CAZy_calls["Firm4_group_1987"]
# 
# 
# 
# panaroo_out[[phylotype]]["Firm4_group_1036", ]
# 
# call <- "GEAANGBI_01530"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "GDAIFKIL_00514"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "HCAMBJLP_01057"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# CAZy_calls_out["Firm4_group_1036", ]
# 
# 
# ####### A few negative ones to confirm:
# 
# "Firm4_whiA" 
# 
# panaroo_out[[phylotype]]["Firm4_whiA", ]
# 
# call <- "PBOALMCJ_00989"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# 
# "Firm4_rplS" 
# 
# panaroo_out[[phylotype]]["Firm4_rplS", ]
# 
# call <- "PBOALMCJ_00730"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# 
# "Firm4_licB"
# panaroo_out[[phylotype]]["Firm4_licB", ]
# 
# call <- "GDAIFKIL_01549"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# 
# # A few other negatives that aren't singletons
# tmp <- panaroo_out[[phylotype]][which(rowSums(panaroo_out[[phylotype]] == "") == 2), ]
# 
# 
# panaroo_out[[phylotype]]["Firm4_group_1573", ]
# #FEIDOIME_00930             NEKIOABH_00775 GEAANGBI_00536              GDAIFKIL_01294              HCAMBJLP_00250
# call <- "FEIDOIME_00930" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "NEKIOABH_00775" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "GEAANGBI_00536" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "GDAIFKIL_01294" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "HCAMBJLP_00250" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# 
# panaroo_out[[phylotype]]["Firm4_group_840", ]
# #FEIDOIME_00920             NEKIOABH_00765              GEAANGBI_00526 GDAIFKIL_01284              HCAMBJLP_00260 
# call <- "FEIDOIME_00920" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "NEKIOABH_00765" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "GEAANGBI_00526" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "GDAIFKIL_01284" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "HCAMBJLP_00260" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# 
# 
# # Check example of 3 annotations to one ortholog.
# #   Bifidobacterium_group_546 GT111,GT2,GT8
# 
# phylotype <- "Bifidobacterium"
# panaroo_out[[phylotype]]["Bifidobacterium_group_546", ]
# # DCJFEAEM_00276 HLODODIH_00264 HLODODIH_00264 HLODODIH_00264
# 
# call <- "DCJFEAEM_00276" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "HLODODIH_00264" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "HLODODIH_00264" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "HLODODIH_00264" 
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# CAZy_calls_out["Bifidobacterium_group_546", ]
# 
# 
# # Some examples with no matching hits?...
# phylotype <- "Snodgrassella"
# panaroo_out[[phylotype]]["Snodgrassella_group_1257", ]
# # DPPCDHFB_00902              HCDKLCFI_01216         NOEGEOHA_01187          FEGGICIM_00330       HFBMLIGJ_01264         
# # BLCANNOI_01324          MEPCKJMN_00968           GICHLMDD_01032        KOOBMECM_00210
# 
# call <- "DPPCDHFB_00902"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "HCDKLCFI_01216"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# call <- "NOEGEOHA_01187"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# 
# call <- "FEGGICIM_00330"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# 
# call <- "HFBMLIGJ_01264"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# 
# call <- "BLCANNOI_01324"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]
# 
# 
# call <- "KOOBMECM_00210"
# call_id <- paste(phylotype, annot2id[[phylotype]][call, "clustering_id"], sep = "_")
# str_split(dbCAN2_all_hits_atleast2[call_id, "Consensus"], ",")[[1]]

