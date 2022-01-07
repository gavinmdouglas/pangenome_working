
rm(list = ls(all.names = TRUE))

draft_id_files <- list.files(path = "/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes_draft", pattern = ".txt", full.names = TRUE)
final_id_files <- list.files(path = "/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes", pattern = ".txt", full.names = TRUE)

draft_ids <- c()
for (f in draft_id_files) {
  draft_ids <- c(draft_ids, read.table(f, stringsAsFactors = FALSE)$V1)
}

final_ids <- c()
for (f in final_id_files) {
  final_ids <- c(final_ids, read.table(f, stringsAsFactors = FALSE)$V1)
}

draft_ids[which(! draft_ids %in% final_ids)]

length(which(! final_ids  %in% draft_ids))

