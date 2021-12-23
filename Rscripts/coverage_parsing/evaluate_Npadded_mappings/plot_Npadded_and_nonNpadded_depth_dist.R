rm(list = ls(all.names = TRUE))

setwd("honey_bee_pangenome/data/Ellegaard_Carrie_pangenome_mapping/Nfilled_mapping_test/")

library(cowplot)
library(ggplot2)

Gilliamella_depth_by_pos_Nfilled <- readRDS(file = "Nfilled_Gilliamella_depth_by_pos.rds")
Gilliamella_depth_by_pos_non.Nfilled <- readRDS(file = "non.Nfilled_Gilliamella_depth_by_pos.rds")

Gilliamella_depth_by_pos_non.Nfilled$xlabels[which(Gilliamella_depth_by_pos_non.Nfilled$xlabels > 0)] <- Gilliamella_depth_by_pos_non.Nfilled$xlabels[which(Gilliamella_depth_by_pos_non.Nfilled$xlabels > 0)] + 50
Gilliamella_depth_by_pos_non.Nfilled$xlabels[which(Gilliamella_depth_by_pos_non.Nfilled$xlabels < 0)] <- Gilliamella_depth_by_pos_non.Nfilled$xlabels[which(Gilliamella_depth_by_pos_non.Nfilled$xlabels < 0)] - 50

leading_dummy <- data.frame(xlabels = 1:50, raw_depth = 0, num_instances = 0)
trailing_dummy <- data.frame(xlabels = -50:-1, raw_depth = 0, num_instances = 0)

Gilliamella_depth_by_pos_non.Nfilled <- rbind(leading_dummy, Gilliamella_depth_by_pos_non.Nfilled)
Gilliamella_depth_by_pos_non.Nfilled <- rbind(Gilliamella_depth_by_pos_non.Nfilled, trailing_dummy)

# Plot ratio of depth between non.Nfilled and Nfilled data:
Gilliamella_depth_by_pos_non.Nfilled_subset <- Gilliamella_depth_by_pos_non.Nfilled
Gilliamella_depth_by_pos_non.Nfilled_subset <- Gilliamella_depth_by_pos_non.Nfilled_subset[-which(Gilliamella_depth_by_pos_non.Nfilled_subset$xlabels %in% c(1:50)), ]
Gilliamella_depth_by_pos_non.Nfilled_subset <- Gilliamella_depth_by_pos_non.Nfilled_subset[-which(Gilliamella_depth_by_pos_non.Nfilled_subset$xlabels %in% c(-50:-1)), ]

Gilliamella_depth_by_pos_Nfilled_subset <- Gilliamella_depth_by_pos_Nfilled
Gilliamella_depth_by_pos_Nfilled_subset <- Gilliamella_depth_by_pos_Nfilled_subset[-which(Gilliamella_depth_by_pos_Nfilled_subset$xlabels %in% c(1:50)), ]
Gilliamella_depth_by_pos_Nfilled_subset <- Gilliamella_depth_by_pos_Nfilled_subset[-which(Gilliamella_depth_by_pos_Nfilled_subset$xlabels %in% c(-50:-1)), ]


Gilliamella_depth_by_pos_non.Nfilled_subset$ratio <- (Gilliamella_depth_by_pos_Nfilled_subset$raw_depth + 1) / (Gilliamella_depth_by_pos_non.Nfilled_subset$raw_depth + 1)
Gilliamella_depth_by_pos_non.Nfilled_subset$log_ratio <- log(Gilliamella_depth_by_pos_non.Nfilled_subset$ratio)

par(mfrow=c(1, 2))
plot(y = Gilliamella_depth_by_pos_non.Nfilled_subset$ratio[which(Gilliamella_depth_by_pos_non.Nfilled_subset$xlabels > 0)],
     x = Gilliamella_depth_by_pos_non.Nfilled_subset$xlabels[which(Gilliamella_depth_by_pos_non.Nfilled_subset$xlabels > 0)],
     ylab = "Depth ratio: (N-filled / non-N-filled)",
     xlab = "Leading gene position")

plot(y = Gilliamella_depth_by_pos_non.Nfilled_subset$ratio[which(Gilliamella_depth_by_pos_non.Nfilled_subset$xlabels < 0)],
     x = Gilliamella_depth_by_pos_non.Nfilled_subset$xlabels[which(Gilliamella_depth_by_pos_non.Nfilled_subset$xlabels < 0)],
     ylab = "Depth ratio: (N-filled / non-N-filled)",
     xlab = "Trailing gene position")

Gilliamella_depth_by_pos_leading_plot_Nfilled <- ggplot(data = Gilliamella_depth_by_pos_Nfilled[which(Gilliamella_depth_by_pos_Nfilled$xlabels > 0), ],
                                                        aes(x = xlabels, y = raw_depth)) +
                                                        geom_bar(stat="identity") +
                                                        theme_bw() +
                                                        geom_vline(xintercept=50, lwd = 0.5, colour = 'red', linetype = "longdash") +
                                                        xlab("Position in gene") +
                                                        ylab("Raw depth") +
                                                        scale_x_continuous(limits=c(0, 301)) +
                                                        ggtitle("N-filled")


Gilliamella_depth_by_pos_leading_plot_non.Nfilled <- ggplot(data = Gilliamella_depth_by_pos_non.Nfilled[which(Gilliamella_depth_by_pos_non.Nfilled$xlabels > 0), ],
                                                            aes(x = xlabels, y = raw_depth)) +
                                                            geom_bar(stat="identity") +
                                                            theme_bw() +
                                                            geom_vline(xintercept=50, lwd = 0.5, colour = 'red', linetype = "longdash") +
                                                            xlab("Position in gene") +
                                                            ylab("Raw depth") +
                                                            scale_x_continuous(limits=c(0, 301)) +
                                                            ggtitle("non-N-filled")


Gilliamella_depth_by_pos_trailing_plot_Nfilled <- ggplot(data = Gilliamella_depth_by_pos_Nfilled[which(Gilliamella_depth_by_pos_Nfilled$xlabels < 0), ],
                                                 aes(x = xlabels, y = raw_depth)) +
                                                geom_bar(stat="identity") +
                                                theme_bw() +
                                                geom_vline(xintercept=-50, lwd = 0.5, colour = 'red', linetype = "longdash") +
                                                xlab("Position in gene") +
                                                ylab("Raw depth") +
                                                ggtitle("N-filled")



Gilliamella_depth_by_pos_trailing_plot_non.Nfilled <- ggplot(data = Gilliamella_depth_by_pos_non.Nfilled[which(Gilliamella_depth_by_pos_non.Nfilled$xlabels < 0), ],
                                                             aes(x = xlabels, y = raw_depth)) +
                                                             geom_bar(stat="identity") +
                                                             theme_bw() +
                                                             geom_vline(xintercept=-50, lwd = 0.5, colour = 'red', linetype = "longdash") +
                                                             xlab("Position in gene") +
                                                             ylab("Raw depth") +
                                                             ggtitle("non-N-filled")


plot_grid(Gilliamella_depth_by_pos_leading_plot_Nfilled,
          Gilliamella_depth_by_pos_leading_plot_non.Nfilled,
          labels = c('a', 'b'), nrow = 2, ncol = 1)


plot_grid(Gilliamella_depth_by_pos_trailing_plot_Nfilled,
          Gilliamella_depth_by_pos_trailing_plot_non.Nfilled,
          labels = c('a', 'b'), nrow = 2, ncol = 1)

