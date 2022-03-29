### Run sanity checks on gene clusters that were produced.

rm(list = ls(all.names = TRUE))

library(pegas)
library(ape)

# obs output
# category        n       S       Wattersons      pi      D
# SRR10810002     374     51      7.845992136447599       13.558931054751904      2.062217025353478
test_fasta = read.dna("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/gene_haplotypes_tajimasD/test_output/SRR10810002_Bifidobacterium_asteroides_yihS.23.seqs.fna",
                      format = "fasta")

nrow(test_fasta)

length(seg.sites(test_fasta))

theta.s(test_fasta)

# Note that the nuc.div function returns pi per site...
# sequences are all 1305 nt
nuc.div(test_fasta) * 1305

tajima.test(test_fasta)$D

