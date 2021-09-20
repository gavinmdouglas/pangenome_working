rm(list = ls(all.names = TRUE))

library("cowplot")
library("ggplot2")

Gilliamella_and_Snodgrassella_breadth <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_and_Snodgrassella_present.rds")
Gilliamella_depth <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Gilliamella_mean_depth_per_site.rds")
Snodgrassella_depth <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/Snodgrassella_mean_depth_per_site.rds")

# Dataset samples
Ellegaard.2019_samples <- read.table("/home/gdouglas/projects/honey_bee/Ellegaard/SRR_ids.txt", stringsAsFactors = FALSE)$V1
Ellegaard.2020_samples <- read.table("/data1/gdouglas/projects/honey_bee/Ellegaard.2020/PRJNA598094_SRR_ids.txt", stringsAsFactors = FALSE)$V1

# Consistently covered core genes
Gilliamella_consistent_core <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/consistently_covered/Gilliamella_2019_and_2020_consistent95_core_panaroo.bed",
                                          header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

Snodgrassella_consistent_core <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/panaroo_breadth_coverage/consistently_covered/Snodgrassella_2019_and_2020_consistent75_core_panaroo.bed",
                                            header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

# Remove all missing rows.
Gilliamella_and_Snodgrassella_breadth <- Gilliamella_and_Snodgrassella_breadth[-which(rowSums(Gilliamella_and_Snodgrassella_breadth) == 0), ]

# Remove rows that are 1 in all samples.
Gilliamella_and_Snodgrassella_breadth <- Gilliamella_and_Snodgrassella_breadth[-which(rowSums(Gilliamella_and_Snodgrassella_breadth == 1) == 74), ]

# Run PCA
Gilliamella_and_Snodgrassella_breadth_PCA <- prcomp(t(Gilliamella_and_Snodgrassella_breadth), scale. = TRUE, center = TRUE)
Gilliamella_and_Snodgrassella_breadth_PCA_df <- data.frame(Gilliamella_and_Snodgrassella_breadth_PCA$x)
Gilliamella_and_Snodgrassella_breadth_PCA_imp <- data.frame(summary(Gilliamella_and_Snodgrassella_breadth_PCA)$importance)
Gilliamella_and_Snodgrassella_breadth_PCA_percent <- data.frame(t(Gilliamella_and_Snodgrassella_breadth_PCA_imp["Proportion of Variance", ])) * 100
Gilliamella_and_Snodgrassella_breadth_PCA_percent$Proportion.of.Variance <- round(Gilliamella_and_Snodgrassella_breadth_PCA_percent$Proportion.of.Variance, digits = 1)
Gilliamella_and_Snodgrassella_breadth_PCA_percent$Proportion.of.Variance <- as.character(Gilliamella_and_Snodgrassella_breadth_PCA_percent$Proportion.of.Variance)

# Country of origin (i.e., which paper it was from)
Gilliamella_and_Snodgrassella_breadth_PCA_df$Country <- NA
Gilliamella_and_Snodgrassella_breadth_PCA_df[Ellegaard.2019_samples, "Country"] <- "Switzerland"
Gilliamella_and_Snodgrassella_breadth_PCA_df[Ellegaard.2020_samples, "Country"] <- "Japan"

# Mean depth per consistently covered core genes (for Gilliamella and Snodgrassella) separately
Gilliamella_and_Snodgrassella_breadth_PCA_df$Gilliamella_core_mean_depth <- colMeans(Gilliamella_depth[Gilliamella_consistent_core, ])[rownames(Gilliamella_and_Snodgrassella_breadth_PCA_df)]
Gilliamella_and_Snodgrassella_breadth_PCA_df$Snodgrassella_core_mean_depth <- colMeans(Snodgrassella_depth[Snodgrassella_consistent_core, ])[rownames(Gilliamella_and_Snodgrassella_breadth_PCA_df)]


# Read in metadata from sample names.
sample_names_2019 <- read.table("/home/gdouglas/projects/honey_bee/Ellegaard/PRJNA473901_SRR_to_sample.txt",
                                header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(sample_names_2019) <- c("SRR", "sample")

sample_names_2019$colony <- gsub("Y.*$", "", sample_names_2019$sample)

sample_names_2019$location <- gsub("_.*$", "", sample_names_2019$sample)
sample_names_2019$location <- gsub("^.*Y", "", sample_names_2019$location)

sample_names_2019$age_raw <- gsub(".*_", "", sample_names_2019$sample)
sample_names_2019$age <- NA
sample_names_2019[grep("W", sample_names_2019$age_raw), "age"] <- "Old"
sample_names_2019[grep("F", sample_names_2019$age_raw), "age"] <- "Middle-aged"
sample_names_2019[grep("N", sample_names_2019$age_raw), "age"] <- "Young"


sample_names_2020 <- read.table("/data1/gdouglas/projects/honey_bee/Ellegaard.2020/PRJNA598094_SRR_to_sample.txt",
                                header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(sample_names_2020) <- c("SRR", "sample")

sample_names_2020$apiary <- NA
sample_names_2020$apiary[grep("Ai", sample_names_2020$sample)] <- "Ai"
sample_names_2020$apiary[grep("Iu", sample_names_2020$sample)] <- "Iu"


Gilliamella_and_Snodgrassella_breadth_PCA_df$Year <- NA
Gilliamella_and_Snodgrassella_breadth_PCA_df[sample_names_2019[which(sample_names_2019$location == "1"), "SRR"], "Year"] <- "2019_1"
Gilliamella_and_Snodgrassella_breadth_PCA_df[sample_names_2019[which(sample_names_2019$location == "2"), "SRR"], "Year"] <- "2019_2"
Gilliamella_and_Snodgrassella_breadth_PCA_df[Ellegaard.2020_samples, "Year"] <- "2020"


Gilliamella_and_Snodgrassella_breadth_PCA_df$Location <- NA
Gilliamella_and_Snodgrassella_breadth_PCA_df[sample_names_2019[which(sample_names_2019$colony == "Dr"), "SRR"], "Location"] <- "Les Droites"
Gilliamella_and_Snodgrassella_breadth_PCA_df[sample_names_2019[which(sample_names_2019$colony == "Gr"), "SRR"], "Location"] <- "Grammont"
Gilliamella_and_Snodgrassella_breadth_PCA_df[sample_names_2020[which(sample_names_2020$apiary == "Ai"), "SRR"], "Location"] <- "Ai"
Gilliamella_and_Snodgrassella_breadth_PCA_df[sample_names_2020[which(sample_names_2020$apiary == "Iu"), "SRR"], "Location"] <- "Iu"

Gilliamella_and_Snodgrassella_breadth_PCA_df$Age <- NA
Gilliamella_and_Snodgrassella_breadth_PCA_df[sample_names_2019$SRR, "Age"] <- sample_names_2019$age
Gilliamella_and_Snodgrassella_breadth_PCA_df[Ellegaard.2020_samples, "Age"] <- "Young"

Gilliamella_and_Snodgrassella_breadth_PCA_df$Age <- factor(Gilliamella_and_Snodgrassella_breadth_PCA_df$Age, levels = c("Young", "Middle-aged", "Old"))

# Age PCA
age_pc1_vs_pc2 <- ggplot(data = Gilliamella_and_Snodgrassella_breadth_PCA_df, aes(x = PC1, y = PC2, shape = Country, colour = Age)) +
                          geom_point(size = 3) +
                          theme_bw() +
                          xlab(paste("PC1 (", Gilliamella_and_Snodgrassella_breadth_PCA_percent["PC1", 1], "%)", sep = "")) +
                          ylab(paste("PC2 (", Gilliamella_and_Snodgrassella_breadth_PCA_percent["PC2", 1], "%)", sep = ""))

age_pc1_vs_pc3 <- ggplot(data = Gilliamella_and_Snodgrassella_breadth_PCA_df, aes(x = PC1, y = PC3, shape = Country, colour = Age)) +
                          geom_point(size = 3) +
                          theme_bw() +
                          xlab(paste("PC1 (", Gilliamella_and_Snodgrassella_breadth_PCA_percent["PC1", 1], "%)", sep = "")) +
                          ylab(paste("PC3 (", Gilliamella_and_Snodgrassella_breadth_PCA_percent["PC3", 1], "%)", sep = ""))

plot_grid(age_pc1_vs_pc2, age_pc1_vs_pc3, nrow = 1, ncol = 2).


# Year PCA
location_pc1_vs_pc2 <- ggplot(data = Gilliamella_and_Snodgrassella_breadth_PCA_df, aes(x = PC1, y = PC2, shape = Country, colour = Year)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab(paste("PC1 (", Gilliamella_and_Snodgrassella_breadth_PCA_percent["PC1", 1], "%)", sep = "")) +
  ylab(paste("PC2 (", Gilliamella_and_Snodgrassella_breadth_PCA_percent["PC2", 1], "%)", sep = ""))

location_pc1_vs_pc3 <- ggplot(data = Gilliamella_and_Snodgrassella_breadth_PCA_df, aes(x = PC1, y = PC3, shape = Country, colour = Year)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab(paste("PC1 (", Gilliamella_and_Snodgrassella_breadth_PCA_percent["PC1", 1], "%)", sep = "")) +
  ylab(paste("PC3 (", Gilliamella_and_Snodgrassella_breadth_PCA_percent["PC3", 1], "%)", sep = ""))

plot_grid(location_pc1_vs_pc2, location_pc1_vs_pc3, nrow = 1, ncol = 2)


# Location PCA
location_pc1_vs_pc2 <- ggplot(data = Gilliamella_and_Snodgrassella_breadth_PCA_df, aes(x = PC1, y = PC2, shape = Country, colour = Location)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab(paste("PC1 (", Gilliamella_and_Snodgrassella_breadth_PCA_percent["PC1", 1], "%)", sep = "")) +
  ylab(paste("PC2 (", Gilliamella_and_Snodgrassella_breadth_PCA_percent["PC2", 1], "%)", sep = ""))

location_pc1_vs_pc3 <- ggplot(data = Gilliamella_and_Snodgrassella_breadth_PCA_df, aes(x = PC1, y = PC3, shape = Country, colour = Location)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab(paste("PC1 (", Gilliamella_and_Snodgrassella_breadth_PCA_percent["PC1", 1], "%)", sep = "")) +
  ylab(paste("PC3 (", Gilliamella_and_Snodgrassella_breadth_PCA_percent["PC3", 1], "%)", sep = ""))

plot_grid(location_pc1_vs_pc2, location_pc1_vs_pc3, nrow = 1, ncol = 2)