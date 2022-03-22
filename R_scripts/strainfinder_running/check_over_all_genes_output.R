# Check over output files to figure out which genes are marked as "passed"
# and which ones have issues to be troubleshooted.

# Also need to figure out which genes have the highest number of tested strains as the best fit: in this case
# additional strain counts should be tried to make sure this is a robust result!

# ONE QUICK FIX WAS THAT I HAD TO RENAME THE Bifidobacterium_coryneforme output folder to be Bifidobacterium_coryneforme_indicum

rm(list = ls(all.names = TRUE))

species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species_present.txt",
                      stringsAsFactors = FALSE)$V1

all_prepped_genes <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/all_genes_prepped.txt",
                                stringsAsFactors = FALSE)$V1

# for (sp in species) {
#   sp_genes <- grep(sp, all_prepped_genes, value = TRUE)
#   
#   print(sp)
#   
#   print(length(sp_genes))
#   
# }

output_dir <- "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_output/"

missing_check <- c()
failed_checks <- list()
max_strains_best_fit <- c()
clean_pass <- c()

num_strains <- c()

for (sp in species) {
  sp_genes <- grep(sp, all_prepped_genes, value = TRUE)
  
  for (gene in sp_genes) {
    
    gene_checkfile <- paste(output_dir, sp, "/", gene, "/", gene, ".check.tsv", sep = "")
    gene_checkfile_op2 <-  paste(output_dir, sp, "/", gene, "/", gene, "_check.tsv", sep = "")
    
    if (file.exists(gene_checkfile)) {
      gene_checkinfo <- read.table(gene_checkfile, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    } else if (file.exists(gene_checkfile_op2)) {
      gene_checkinfo <- read.table(gene_checkfile_op2, sep = "\t", header = FALSE, stringsAsFactors = FALSE) 
    } else {
      missing_check <- c(missing_check, gene)
      next
    }
    
    if (gene_checkinfo$V2 != "passed") {
      print(gene)
     if (gene_checkinfo$V2 %in% names(failed_checks)) {
       failed_checks[[gene_checkinfo$V2]] <- c(failed_checks[[gene_checkinfo$V2]], gene_checkinfo$V2)
     } else {
       failed_checks[[gene_checkinfo$V2]] <- gene_checkinfo$V2
     }
    } else if(gene_checkinfo$V3 >= 29) {
      max_strains_best_fit <- c(max_strains_best_fit, gene)
      num_strains <- c(num_strains, gene_checkinfo$V3)
    } else {
      clean_pass <- c(clean_pass, gene)
      num_strains <- c(num_strains, gene_checkinfo$V3)
    }
    
  }
}

length(clean_pass)
length(max_strains_best_fit)
length(missing_check)
failed_checks

write.table(x = max_strains_best_fit, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/working_gene_sets/batch3_29_or_30_best_fit.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)




# Check second batch of ids.

rm(list = ls(all.names = TRUE))

species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species.txt",
                      stringsAsFactors = FALSE)$V1

all_prepped_genes <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/working_gene_sets/batch3_29_or_30_best_fit.txt",
                                stringsAsFactors = FALSE)$V1

output_dir <- "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_output_more_strains/"

missing_check <- c()
failed_checks <- list()
max_strains_best_fit <- c()
min_strains_best_fit <- c()
clean_pass <- c()

num_strains <- c()

for (sp in species) {
  sp_genes <- grep(sp, all_prepped_genes, value = TRUE)
  
  for (gene in sp_genes) {
    
    gene_checkfile <- paste(output_dir, sp, "/", gene, "/", gene, ".check.tsv", sep = "")
    gene_checkfile_op2 <-  paste(output_dir, sp, "/", gene, "/", gene, "_check.tsv", sep = "")
    
    if (file.exists(gene_checkfile)) {
      gene_checkinfo <- read.table(gene_checkfile, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    } else if (file.exists(gene_checkfile_op2)) {
      gene_checkinfo <- read.table(gene_checkfile_op2, sep = "\t", header = FALSE, stringsAsFactors = FALSE) 
    } else {
      missing_check <- c(missing_check, gene)
      next
    }
    
    if (gene_checkinfo$V2 != "passed") {
      print(gene)
      if (gene_checkinfo$V2 %in% names(failed_checks)) {
        failed_checks[[gene_checkinfo$V2]] <- c(failed_checks[[gene_checkinfo$V2]], gene_checkinfo$V2)
      } else {
        failed_checks[[gene_checkinfo$V2]] <- gene_checkinfo$V2
      }
    } else if(gene_checkinfo$V3 >= 73) {
      max_strains_best_fit <- c(max_strains_best_fit, gene)
      num_strains <- c(num_strains, gene_checkinfo$V3)
    } else if(gene_checkinfo$V3 <= 30) {
      min_strains_best_fit <- c(min_strains_best_fit, gene)
      num_strains <- c(num_strains, gene_checkinfo$V3)  
    
    } else {
      clean_pass <- c(clean_pass, gene)
      num_strains <- c(num_strains, gene_checkinfo$V3)
    }
    
  }
}

length(clean_pass)
length(max_strains_best_fit)
length(missing_check)
failed_checks

write.table(x = clean_pass, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/all_gene_input/working_gene_sets/batch3_genes_to_copy.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)
