### Heatmap of species presence/absence with sample metadata information shown as well.

rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(ggplot2)

species_25per <- read.table(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/summary/species_presence_core_25percent.tsv.gz",
                            row.names = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

species_90per <- read.table(file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/summary/species_presence_core_90percent.tsv.gz",
                            row.names = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

sample_metadata <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard.2019.2020_metadata.tsv",
                              header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

column_annot <- HeatmapAnnotation(Country=sample_metadata$Country,
                                  Apiary=sample_metadata$Apiary,
                                  Year=sample_metadata$Year,
                                  Age=factor(sample_metadata$Age, levels = c("Young", "Middle-aged", "Old")))

species_25per <- t(species_25per[rownames(sample_metadata), ])
species_90per <- t(species_90per[rownames(sample_metadata), ])

Heatmap(as.matrix(species_25per),
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        clustering_distance_rows = "binary",
        clustering_distance_columns = "binary",
        top_annotation = column_annot)


Heatmap(as.matrix(species_90per),
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        clustering_distance_rows = "binary",
        clustering_distance_columns = "binary",
        top_annotation = column_annot)

