### The script Validate.pl that is part of the CONCOCT directory expects bins to be numeric names
### Needed to rename the ids for the das tool, maxbin2, and metabat2 outputs.

rm(list = ls(all.names = TRUE))

setwd("/Users/Gavin/Google_Drive/postdoc/honey_bee_pangenome/data/binning_output/")

dastool_scaffold2bin <- read.table("megahit_cosassembly_scaffold2bin_files/output_DASTool_scaffolds2bin.txt",
                                   header = FALSE, sep = "\t", stringsAsFactors = FALSE)
dastool_scaffold2bin$V2 <- as.numeric(factor(dastool_scaffold2bin$V2))
write.table(x = dastool_scaffold2bin,
            file = "concoct_validate_prepped_files/output_DASTool_scaffolds2bin_as_numeric.csv",
            sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

maxbin2_scaffold2bin <- read.table("megahit_cosassembly_scaffold2bin_files/maxbin2_scaffold2bin.tsv",
                                   header = FALSE, sep = "\t", stringsAsFactors = FALSE)
maxbin2_scaffold2bin$V2 <- as.numeric(factor(maxbin2_scaffold2bin$V2))
write.table(x = maxbin2_scaffold2bin,
            file = "concoct_validate_prepped_files/output_maxbin2_scaffold2bin_as_numeric.csv",
            sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)


metabat2_scaffold2bin <- read.table("megahit_cosassembly_scaffold2bin_files/metabat2_scaffold2bin.tsv",
                                   header = FALSE, sep = "\t", stringsAsFactors = FALSE)
metabat2_scaffold2bin$V2 <- as.numeric(factor(metabat2_scaffold2bin$V2))
write.table(x = metabat2_scaffold2bin,
            file = "concoct_validate_prepped_files/metabat2_scaffold2bin_as_numeric.csv",
            sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

