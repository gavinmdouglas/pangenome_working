### Check over functional annotation files

rm(list = ls(all.names = TRUE))

library(stringr)

### 1. Check over a few lines of the clean master files with raw files manually

# Read in master func files for mapping to CDS.
eggNOG_output <- read.table(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/all_microbiota_eggNOG.tsv",
                            header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, quote = "", comment.char = "")

KO_output <- read.table(file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/KEGG_output/all_microbiota_KO_calls.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)

# These lines are from the eggNOG raw files (e.g., all_microbiota_proteins.part-2_eggNOG_mapper_out.txt.gz)
# BKLEFNHF_00756 1343158.SACS_0112       2.96E-124       358     "COG0491@1|root,COG0491@2|Bacteria,1MU8Q@1224|Proteobacteria,2TSVS@28211|Alphaproteobacteria"
# IEOHIMPA_01448  1545701.LACWKB10_0131   0       1261    "COG0542@1|root,COG0542@2|Bacteria,1TPMU@1239|Firmicutes,4HA0V@91061|Bacilli,3F3K9@33958|Lactobacillaceae"      91061|Bacilli   O       Belongs to the ClpA ClpB family clpE    -       -       ko:K03697       -       -       -       -       "ko00000,ko03110"       -       -       -       "AAA,AAA_2,ClpB_D2-small,UVR"
# BKLEFNHF_00707  663932.KB902575_gene1113        1.51E-49        166     "COG5342@1|root,COG5342@2|Bacteria"     2|Bacteria      S       invasion associated locus B     -       -       -       -       -       -       -       -       -       -       -       -       IalB
eggNOG_output[c("BKLEFNHF_00756", "IEOHIMPA_01448", "BKLEFNHF_00707"), ]


# From KO raw file (user_ko.txt.gz)
# IDCHIMDP_02173
# IDCHIMDP_02174  K00549
# HHGPBACN_00777  K17077
# HHGPBACN_00777  K23059
is.na(KO_output["IDCHIMDP_02173", ])
KO_output["IDCHIMDP_02174", ] == "K00549"
KO_output["HHGPBACN_00777", ] == "K17077,K23059"



### 2. Check on some CAZy's by hand.
### Pay special attention to multi-annotations as these originally had quotes around them in the output.
### Check these lines of the output line.
# Bifidobacterium_asteroides_group_706    GH43,GH51
# Apilactobacillus_kunkeei_group_847      GH13
# Gilliamella_apicola_group_2271  GT2,GT4
# Also confirm that a couple are NOT present in CAZy output file makes sense.

CAZy_output <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_output/all_microbiota_CAZy_by_ortholog.tsv",
                          header = FALSE, row.names = 1, stringsAsFactors = FALSE)

# Bifidobacterium_asteroides_group_706    GH43,GH51
Bifidobacterium_asteroides_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Bifidobacterium_asteroides/gene_presence_absence_roary.csv.gz",
                                                 header = TRUE, stringsAsFactors = FALSE, sep = ",", comment.char = "", quote = "")
rownames(Bifidobacterium_asteroides_panaroo) <- paste("Bifidobacterium_asteroides", Bifidobacterium_asteroides_panaroo$Gene, sep = "_")
Bifidobacterium_asteroides_test_cds_ids <- as.character(Bifidobacterium_asteroides_panaroo["Bifidobacterium_asteroides_group_706", 15:ncol(Bifidobacterium_asteroides_panaroo)])
if (length(which(Bifidobacterium_asteroides_test_cds_ids == "")) > 0) {
  Bifidobacterium_asteroides_test_cds_ids <- Bifidobacterium_asteroides_test_cds_ids[-which(Bifidobacterium_asteroides_test_cds_ids == "")]
}
eggNOG_output[Bifidobacterium_asteroides_test_cds_ids, "CAZy"]


# Apilactobacillus_kunkeei_group_847      GH13
Apilactobacillus_kunkeei_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Apilactobacillus_kunkeei/gene_presence_absence_roary.csv.gz",
                                                 header = TRUE, stringsAsFactors = FALSE, sep = ",", comment.char = "", quote = "")
rownames(Apilactobacillus_kunkeei_panaroo) <- paste("Apilactobacillus_kunkeei", Apilactobacillus_kunkeei_panaroo$Gene, sep = "_")
Apilactobacillus_kunkeei_test_cds_ids <- as.character(Apilactobacillus_kunkeei_panaroo["Apilactobacillus_kunkeei_group_847", 15:ncol(Apilactobacillus_kunkeei_panaroo)])
if (length(which(Apilactobacillus_kunkeei_test_cds_ids == "")) > 0) {
  Apilactobacillus_kunkeei_test_cds_ids <- Apilactobacillus_kunkeei_test_cds_ids[-which(Apilactobacillus_kunkeei_test_cds_ids == "")]
}
eggNOG_output[Apilactobacillus_kunkeei_test_cds_ids, "CAZy"]



# Gilliamella_apicola_group_2271  GT2,GT4
Gilliamella_apicola_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Gilliamella_apicola/gene_presence_absence_roary.csv.gz",
                                               header = TRUE, stringsAsFactors = FALSE, sep = ",", comment.char = "", quote = "")
rownames(Gilliamella_apicola_panaroo) <- paste("Gilliamella_apicola", Gilliamella_apicola_panaroo$Gene, sep = "_")
Gilliamella_apicola_test_cds_ids <- as.character(Gilliamella_apicola_panaroo["Gilliamella_apicola_group_2271", 15:ncol(Gilliamella_apicola_panaroo)])
if (length(which(Gilliamella_apicola_test_cds_ids == "")) > 0) {
  Gilliamella_apicola_test_cds_ids <- Gilliamella_apicola_test_cds_ids[-which(Gilliamella_apicola_test_cds_ids == "")]
}
eggNOG_output[Gilliamella_apicola_test_cds_ids, "CAZy"]


# Several not found in CAZy
missing_CAZy_test_ids1 <- as.character(Gilliamella_apicola_panaroo["Gilliamella_apicola_flgA", 15:ncol(Gilliamella_apicola_panaroo)])
missing_CAZy_test_ids2 <- as.character(Gilliamella_apicola_panaroo["Gilliamella_apicola_group_2543", 15:ncol(Gilliamella_apicola_panaroo)])
missing_CAZy_test_ids3 <- as.character(Gilliamella_apicola_panaroo["Gilliamella_apicola_cdd", 15:ncol(Gilliamella_apicola_panaroo)])

eggNOG_output[missing_CAZy_test_ids1, "CAZy"]
eggNOG_output[missing_CAZy_test_ids2, "CAZy"]
eggNOG_output[missing_CAZy_test_ids3, "CAZy"]




### 3. Check on some COGs by hand.
### Check these lines of the output line.
# Bifidobacterium_asteroides_group_706    G
# Apilactobacillus_apinorum_DPLMPBNN_00008        E,G,P
# Gilliamella_apicola_barS1       I,Q

# Bifidobacterium_asteroides_group_706    G
Bifidobacterium_asteroides_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Bifidobacterium_asteroides/gene_presence_absence_roary.csv.gz",
                                                 header = TRUE, stringsAsFactors = FALSE, sep = ",", comment.char = "", quote = "")
rownames(Bifidobacterium_asteroides_panaroo) <- paste("Bifidobacterium_asteroides", Bifidobacterium_asteroides_panaroo$Gene, sep = "_")
Bifidobacterium_asteroides_test_cds_ids <- as.character(Bifidobacterium_asteroides_panaroo["Bifidobacterium_asteroides_group_706", 15:ncol(Bifidobacterium_asteroides_panaroo)])
if (length(which(Bifidobacterium_asteroides_test_cds_ids == "")) > 0) {
  Bifidobacterium_asteroides_test_cds_ids <- Bifidobacterium_asteroides_test_cds_ids[-which(Bifidobacterium_asteroides_test_cds_ids == "")]
}
eggNOG_output[Bifidobacterium_asteroides_test_cds_ids, "COG_category"]


# Apilactobacillus_apinorum_DPLMPBNN_00008        E,G,P
Apilactobacillus_apinorum_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Apilactobacillus_apinorum/gene_presence_absence_roary.csv.gz",
                                                 header = TRUE, stringsAsFactors = FALSE, sep = ",", comment.char = "", quote = "")
rownames(Apilactobacillus_apinorum_panaroo) <- paste("Apilactobacillus_apinorum", Apilactobacillus_apinorum_panaroo$Gene, sep = "_")
Apilactobacillus_apinorum_test_cds_ids <- as.character(Apilactobacillus_apinorum_panaroo["Apilactobacillus_apinorum_DPLMPBNN_00008", 15:ncol(Apilactobacillus_apinorum_panaroo)])
if (length(which(Apilactobacillus_apinorum_test_cds_ids == "")) > 0) {
  Apilactobacillus_apinorum_test_cds_ids <- Apilactobacillus_apinorum_test_cds_ids[-which(Apilactobacillus_apinorum_test_cds_ids == "")]
}
eggNOG_output[Apilactobacillus_apinorum_test_cds_ids, "COG_category"]


# Gilliamella_apicola_barS1       I,Q
Gilliamella_apicola_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Gilliamella_apicola/gene_presence_absence_roary.csv.gz",
                                          header = TRUE, stringsAsFactors = FALSE, sep = ",", comment.char = "", quote = "")
rownames(Gilliamella_apicola_panaroo) <- paste("Gilliamella_apicola", Gilliamella_apicola_panaroo$Gene, sep = "_")
Gilliamella_apicola_test_cds_ids <- as.character(Gilliamella_apicola_panaroo["Gilliamella_apicola_barS1", 15:ncol(Gilliamella_apicola_panaroo)])
if (length(which(Gilliamella_apicola_test_cds_ids == "")) > 0) {
  Gilliamella_apicola_test_cds_ids <- Gilliamella_apicola_test_cds_ids[-which(Gilliamella_apicola_test_cds_ids == "")]
}
eggNOG_output[Gilliamella_apicola_test_cds_ids, "COG_category"]




### 4. Check on some broadest OGs by hand.
### Check these lines of the output line.
# Bifidobacterium_asteroides_group_706    COG3507
# Apilactobacillus_apinorum_DPLMPBNN_00008        COG0477
# Gilliamella_apicola_barS1       COG1028


# Bifidobacterium_asteroides_group_706    COG3507
Bifidobacterium_asteroides_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Bifidobacterium_asteroides/gene_presence_absence_roary.csv.gz",
                                                 header = TRUE, stringsAsFactors = FALSE, sep = ",", comment.char = "", quote = "")
rownames(Bifidobacterium_asteroides_panaroo) <- paste("Bifidobacterium_asteroides", Bifidobacterium_asteroides_panaroo$Gene, sep = "_")
Bifidobacterium_asteroides_test_cds_ids <- as.character(Bifidobacterium_asteroides_panaroo["Bifidobacterium_asteroides_group_706", 15:ncol(Bifidobacterium_asteroides_panaroo)])
if (length(which(Bifidobacterium_asteroides_test_cds_ids == "")) > 0) {
  Bifidobacterium_asteroides_test_cds_ids <- Bifidobacterium_asteroides_test_cds_ids[-which(Bifidobacterium_asteroides_test_cds_ids == "")]
}
eggNOG_output[Bifidobacterium_asteroides_test_cds_ids, "broadest_OG"]


# Apilactobacillus_apinorum_DPLMPBNN_00008        COG0477
Apilactobacillus_apinorum_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Apilactobacillus_apinorum/gene_presence_absence_roary.csv.gz",
                                                header = TRUE, stringsAsFactors = FALSE, sep = ",", comment.char = "", quote = "")
rownames(Apilactobacillus_apinorum_panaroo) <- paste("Apilactobacillus_apinorum", Apilactobacillus_apinorum_panaroo$Gene, sep = "_")
Apilactobacillus_apinorum_test_cds_ids <- as.character(Apilactobacillus_apinorum_panaroo["Apilactobacillus_apinorum_DPLMPBNN_00008", 15:ncol(Apilactobacillus_apinorum_panaroo)])
if (length(which(Apilactobacillus_apinorum_test_cds_ids == "")) > 0) {
  Apilactobacillus_apinorum_test_cds_ids <- Apilactobacillus_apinorum_test_cds_ids[-which(Apilactobacillus_apinorum_test_cds_ids == "")]
}
eggNOG_output[Apilactobacillus_apinorum_test_cds_ids, "broadest_OG"]


# Gilliamella_apicola_barS1       COG1028
Gilliamella_apicola_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Gilliamella_apicola/gene_presence_absence_roary.csv.gz",
                                          header = TRUE, stringsAsFactors = FALSE, sep = ",", comment.char = "", quote = "")
rownames(Gilliamella_apicola_panaroo) <- paste("Gilliamella_apicola", Gilliamella_apicola_panaroo$Gene, sep = "_")
Gilliamella_apicola_test_cds_ids <- as.character(Gilliamella_apicola_panaroo["Gilliamella_apicola_barS1", 15:ncol(Gilliamella_apicola_panaroo)])
if (length(which(Gilliamella_apicola_test_cds_ids == "")) > 0) {
  Gilliamella_apicola_test_cds_ids <- Gilliamella_apicola_test_cds_ids[-which(Gilliamella_apicola_test_cds_ids == "")]
}
eggNOG_output[Gilliamella_apicola_test_cds_ids, "broadest_OG"]





### 4. Check on some KOs by hand.
### Check these lines of the output line.
# Bifidobacterium_asteroides_group_706    K01198
# Apilactobacillus_kunkeei_macB   K02003,K02004

# Also confirm that a couple are NOT present in output file makes sense.
## Apilactobacillus_apinorum_DPLMPBNN_00008        --> missing
## Gilliamella_apicola_barS1 --> missing

# Bifidobacterium_asteroides_group_706    K01198
Bifidobacterium_asteroides_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Bifidobacterium_asteroides/gene_presence_absence_roary.csv.gz",
                                                 header = TRUE, stringsAsFactors = FALSE, sep = ",", comment.char = "", quote = "")
rownames(Bifidobacterium_asteroides_panaroo) <- paste("Bifidobacterium_asteroides", Bifidobacterium_asteroides_panaroo$Gene, sep = "_")
Bifidobacterium_asteroides_test_cds_ids <- as.character(Bifidobacterium_asteroides_panaroo["Bifidobacterium_asteroides_group_706", 15:ncol(Bifidobacterium_asteroides_panaroo)])
if (length(which(Bifidobacterium_asteroides_test_cds_ids == "")) > 0) {
  Bifidobacterium_asteroides_test_cds_ids <- Bifidobacterium_asteroides_test_cds_ids[-which(Bifidobacterium_asteroides_test_cds_ids == "")]
}
KO_output[Bifidobacterium_asteroides_test_cds_ids, 1]


# Apilactobacillus_kunkeei_macB   K02003,K02004
Apilactobacillus_kunkeei_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Apilactobacillus_kunkeei/gene_presence_absence_roary.csv.gz",
                                                header = TRUE, stringsAsFactors = FALSE, sep = ",", comment.char = "", quote = "")
rownames(Apilactobacillus_kunkeei_panaroo) <- paste("Apilactobacillus_kunkeei", Apilactobacillus_kunkeei_panaroo$Gene, sep = "_")
Apilactobacillus_kunkeei_test_cds_ids <- as.character(Apilactobacillus_kunkeei_panaroo["Apilactobacillus_kunkeei_macB", 15:ncol(Apilactobacillus_kunkeei_panaroo)])
if (length(which(Apilactobacillus_kunkeei_test_cds_ids == "")) > 0) {
  Apilactobacillus_kunkeei_test_cds_ids <- Apilactobacillus_kunkeei_test_cds_ids[-which(Apilactobacillus_kunkeei_test_cds_ids == "")]
}
KO_output[Apilactobacillus_kunkeei_test_cds_ids, 1]


## Apilactobacillus_apinorum_DPLMPBNN_00008        --> missing
## Gilliamella_apicola_barS1 --> missing

Apilactobacillus_apinorum_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Apilactobacillus_apinorum/gene_presence_absence_roary.csv.gz",
                                               header = TRUE, stringsAsFactors = FALSE, sep = ",", comment.char = "", quote = "")
rownames(Apilactobacillus_apinorum_panaroo) <- paste("Apilactobacillus_apinorum", Apilactobacillus_apinorum_panaroo$Gene, sep = "_")
Apilactobacillus_apinorum_test_cds_ids <- as.character(Apilactobacillus_apinorum_panaroo["Apilactobacillus_apinorum_DPLMPBNN_00008", 15:ncol(Apilactobacillus_apinorum_panaroo)])
if (length(which(Apilactobacillus_apinorum_test_cds_ids == "")) > 0) {
  Apilactobacillus_apinorum_test_cds_ids <- Apilactobacillus_apinorum_test_cds_ids[-which(Apilactobacillus_apinorum_test_cds_ids == "")]
}
KO_output[Apilactobacillus_apinorum_test_cds_ids, 1]



Gilliamella_apicola_panaroo <- read.table("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/Gilliamella_apicola/gene_presence_absence_roary.csv.gz",
                                                header = TRUE, stringsAsFactors = FALSE, sep = ",", comment.char = "", quote = "")
rownames(Gilliamella_apicola_panaroo) <- paste("Gilliamella_apicola", Gilliamella_apicola_panaroo$Gene, sep = "_")
Gilliamella_apicola_test_cds_ids <- as.character(Gilliamella_apicola_panaroo["Gilliamella_apicola_barS1", 15:ncol(Gilliamella_apicola_panaroo)])
if (length(which(Gilliamella_apicola_test_cds_ids == "")) > 0) {
  Gilliamella_apicola_test_cds_ids <- Gilliamella_apicola_test_cds_ids[-which(Gilliamella_apicola_test_cds_ids == "")]
}
KO_output[Gilliamella_apicola_test_cds_ids, 1]