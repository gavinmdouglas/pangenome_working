rm(list = ls(all.names = TRUE))

# Check that the mean coverage reported by bedtools coverage is indeed the mean of all per-site depths.

setwd("honey_bee_pangenome/data/Ellegaard_panaroo_pangenome_mapping/tests/test_mean_gene_coverage")

mean_out <- read.table("SRR7287233.Gilliamella.merged.nonparalog.mean_genes.TEST.bedGraph.gz",
                       header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

per.site_out <- read.table("SRR7287233.Gilliamella.merged.nonparalog.per_site.bedGraph.gz",
                            header = FALSE, sep = "\t")

per.site_out <- per.site_out[, c(-2, -3, -4)]

per.site_out_mean <- aggregate(V5 ~ V1, data = per.site_out, FUN = mean)

rownames(per.site_out_mean) <- per.site_out_mean$V1
per.site_out_mean <- per.site_out_mean[rownames(mean_out), , drop = FALSE]

ll.equal(per.site_out_mean$V5, mean_out$V4)
###[1] "Mean relative difference: 2.121736e-08" (consistent with rounding differences).

