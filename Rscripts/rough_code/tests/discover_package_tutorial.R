library(discover)

data("BRCA.mut")

BRCA.mut[1:5, 1:5]

events <- discover.matrix(BRCA.mut, strata = sample(c("group1", "group2"), size = ncol(BRCA.mut), replace = TRUE))

subset <- rowSums(BRCA.mut) > 25

result.mutex <- pairwise.discover.test(events[subset, ], fdr.method = "DBH", alternative = "less")
result.mutex_df <- as.data.frame(result.mutex, q.threshold = 0.01)

result.cooccur <- pairwise.discover.test(events[subset, ], fdr.method = "DBH", alternative = "greater")
result.cooccur_df <- as.data.frame(result.cooccur, q.threshold = 0.05)
