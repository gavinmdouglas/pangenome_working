rm(list = ls(all.names = TRUE))

# As many of these genomes were already downloaded as part of the original batch, need to figure out what accessions need to be downloaded.
# Need to also figure out how many accessions in original batch don't overlap to get an idea of how different the genome sets really are.

generate_rsync_download_cmds <- function(in_links, outdir) {
  
  command_start = "rsync --copy-links --times "
  fasta_suffix <- "_genomic.fna.gz"
  
  cmds <- c()
  
  for(l in in_links) {
    l <- sub("^ftp", "rsync", l)
    link_basename <- basename(l)
    cmds <- c(cmds, paste0(command_start, l, "/", link_basename, fasta_suffix, " ", outdir, "/", link_basename, fasta_suffix))
    
  }
  
  return(cmds)
}

generate_rsync_download_cmds_md5sum <- function(in_links, outdir) {
  
  command_start = "rsync --copy-links --times "
  fasta_suffix <- "_genomic.fna.gz"
  
  cmds <- c()
  
  for(l in in_links) {
    l <- sub("^ftp", "rsync", l)
    link_basename <- basename(l)
    cmds <- c(cmds, paste0(command_start, l, "/md5checksums.txt ", outdir, "/", link_basename, fasta_suffix, ".exp.md5"))
    
  }
  
  return(cmds)
}

phylotypes <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/phylotypes.txt", stringsAsFactors = FALSE)$V1


phylotype_orig_ids <- list()
orig_ids <- c()

for (p in phylotypes) {
  phylotype_id_filename <- paste0("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_05_20_Carrie_genome_accessions_renamed/", p, "_ids.txt")
  phylotype_orig_ids[[p]] <- read.table(phylotype_id_filename, stringsAsFactors = FALSE)$V1
  phylotype_orig_ids[[p]] <- gsub("^GCA_", "GCA-", phylotype_orig_ids[[p]])
  phylotype_orig_ids[[p]] <- gsub("^GCF_", "GCF-", phylotype_orig_ids[[p]])
  phylotype_orig_ids[[p]] <- gsub("_.*$", "", phylotype_orig_ids[[p]])
  phylotype_orig_ids[[p]] <- gsub("^GCA-", "GCA_", phylotype_orig_ids[[p]])
  phylotype_orig_ids[[p]] <- gsub("^GCF-", "GCF_", phylotype_orig_ids[[p]])
  
  orig_ids <- c(orig_ids, phylotype_orig_ids[[p]])
}

# Remove everything about genome accession (used quick hack approach that required changing the prefix first)


# Read in all accessions that want to analyze.

species_names <- gsub(".tsv$", "", list.files("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process/"))
accession_files <- list.files("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/accessions_to_process", full.names = TRUE)

compare_log <- data.frame(matrix(NA, nrow = length(species_names), ncol = 2))
rownames(compare_log) <- species_names
colnames(compare_log) <- c("total", "download_needed")

mkdir_cmds <- paste("mkdir", species_names)

download_cmds <- c()
md5_download_cmds <- c()

all_new_accessions <- c()

for (i in 1:length(accession_files)) {

  accession_file <- accession_files[i]
  species_name <- species_names[i]
  
  info_tab <- read.table(accession_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  all_new_accessions <- c(all_new_accessions, info_tab$accession)
  
  accessions_i_download_needed <- which(! info_tab$accession %in% orig_ids)
  
  compare_log[species_name, ] <- c(nrow(info_tab), length(accessions_i_download_needed))
  
  info_tab_to_download <- info_tab[accessions_i_download_needed, ]
  
  download_cmds <- c(download_cmds, generate_rsync_download_cmds(info_tab_to_download$NCBI_download, species_name))

  md5_download_cmds <- c(md5_download_cmds, generate_rsync_download_cmds_md5sum(info_tab_to_download$NCBI_download, species_name))
}

write.table(x = mkdir_cmds, file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/new_genome_downloads/mkdir_cmds.sh",
            sep = "", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(x = download_cmds, file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/new_genome_downloads/download_cmds.sh",
            sep = "", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(x = md5_download_cmds, file = "/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/new_genome_downloads/md5_download_cmds.sh",
            sep = "", col.names = FALSE, row.names = FALSE, quote = FALSE)


# Also figure out how many genomes you already have will not be used.

for (p in phylotypes) {
 print(p)
  print(length(phylotype_orig_ids[[p]]))
  print(length(which(! phylotype_orig_ids[[p]] %in% all_new_accessions)))
}


### Used these commands to figure out which accessions were in the original set and were missed when I did a custom scan

bifido <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_05_20_Carrie_genome_accessions/bifido_strains_host.txt",
                     header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

bifido <- bifido[grep("mellifera", bifido)]
bifido <- gsub("^GCA_", "GCA-", bifido)
bifido <- gsub("^GCF_", "GCF-", bifido)
bifido <- gsub("_.*$", "", bifido)
bifido <- gsub("^GCA-", "GCA_", bifido)
bifido <- gsub("^GCF-", "GCF_", bifido)
bifido[which(! bifido %in% all_new_accessions)]


# asteroides to add
# GCA_002715865.1
# GCA_002846895.1
# "GCA_007559155.1"
# "GCA_007559275.1"




firm4 <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_05_20_Carrie_genome_accessions/firm4_strains_host.txt",
                     header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

firm4 <- firm4[grep("mellifera", firm4)]
firm4 <- gsub("^GCA_", "GCA-", firm4)
firm4 <- gsub("^GCF_", "GCF-", firm4)
firm4 <- gsub("_.*$", "", firm4)
firm4 <- gsub("^GCA-", "GCA_", firm4)
firm4 <- gsub("^GCF-", "GCF_", firm4)
firm4[which(! firm4 %in% all_new_accessions)]

# All unusual Lactobacillus species that didn't intersect with the species I parsed, but I guess are in Firm4


firm5 <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_05_20_Carrie_genome_accessions/lacto5_strains_host.txt",
                    header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

firm5 <- firm5[grep("mellifera", firm5)]
firm5 <- gsub("^GCA_", "GCA-", firm5)
firm5 <- gsub("^GCF_", "GCF-", firm5)
firm5 <- gsub("_.*$", "", firm5)
firm5 <- gsub("^GCA-", "GCA_", firm5)
firm5 <- gsub("^GCF-", "GCF_", firm5)
firm5[which(! firm5 %in% all_new_accessions)]

# Lactobacillus apis
# GCA_002837055.1
# GCA_900094785.1

gilli <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_05_20_Carrie_genome_accessions/gilli_strains_host.txt",
                    header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

gilli <- gilli[grep("mellifera", gilli)]
gilli <- gsub("^GCA_", "GCA-", gilli)
gilli <- gsub("^GCF_", "GCF-", gilli)
gilli <- gsub("_.*$", "", gilli)
gilli <- gsub("^GCA-", "GCA_", gilli)
gilli <- gsub("^GCF-", "GCF_", gilli)
gilli[which(! gilli %in% all_new_accessions)]

# Gilliamella apicola
# GCA_000599985.1


snod <- read.table("/data1/gdouglas/projects/honey_bee/adding_new_microbiota_genomes/2021_05_20_Carrie_genome_accessions/snod_strains_host.txt",
                   header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

snod <- snod[grep("mellifera", snod)]
snod <- gsub("^GCA_", "GCA-", snod)
snod <- gsub("^GCF_", "GCF-", snod)
snod <- gsub("_.*$", "", snod)
snod <- gsub("^GCA-", "GCA_", snod)
snod <- gsub("^GCF-", "GCF_", snod)
snod[which(! snod %in% all_new_accessions)]



