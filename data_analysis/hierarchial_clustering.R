library(cluster)
library(gplots)
library(ggdendro)
library(ggplot2)
library(grid)
library(gridExtra)
library(tidyverse)

setwd("/projectnb/bf528/users/hedgehog/project_1/noise_filtering/")

filtered <- read.csv("hedgehog1_data-filter_all.csv")
filtered <- filtered[,-1] # remove first column of garbage

# turn file names into probe set names
raw_names <- strsplit(names(filtered),'_')
names(filtered) <- apply(data.frame(names(filtered)), 1, 
                         function(x) strsplit(x, "_")[[1]][1])

# 5.1 + 5.2
# find distance between points and then cut tree
d <- dist(t(filtered[,-1]), method = "euclidean")
fit <- hclust(d, method = "ward.D")
# check agglomerative coeffecient, higher is better
agnes_result <- agnes(t(filtered[,-1]), method="average")$ac
print(paste0("Agnes agg. coefficient: ", agnes_result, ", higher is better."))
groups <- cutree(fit, k = 2)  # cut tree into 2 clusters
print(table(groups)) # number of samples in each group
fit <- as.dendrogram(fit)

# put way too much work into plotting:
upper <- ggdendrogram(fit, theme_dendro = FALSE, cex = 0) + 
  xlab("Samples") + 
  ylab("Height") +
  ggtitle("a)") +
  annotate("rect", xmin = 0, xmax = 77.5, ymin = 0, ymax = 210, 
           fill="red", colour="red", alpha=0.1) +
  # annotate("text", x = 3, y = 250, label = "b)") +
  annotate("rect", xmin = 77.6, xmax = 134, ymin = 0, ymax = 210, 
           fill="orange", colour="orange", alpha=0.1) +
  # annotate("text", x = 80, y = 250, label = "c)") +
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="black" ),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

lower_1 <- ggdendrogram(cut(fit, 500)$lower[[1]], theme_dendro = FALSE) + 
  xlab("Sample name") + 
  ylab("Height") +
  ylim(c(0, 210)) +
  ggtitle("b)") +
  theme(panel.grid.major.y = element_line(size=.1, color="black", 
                                          linetype = "dashed"),
        panel.background = element_rect(fill = alpha("red", 0.1),
                                        colour = alpha("red", 0.1),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
lower_2 <- ggdendrogram(cut(fit, 500)$lower[[2]], theme_dendro = FALSE) + 
  xlab("Sample name") + 
  ylab("Height") +
  ylim(c(0, 210)) +
  ggtitle("c)") +
  theme(panel.grid.major.y = element_line(size=.1, color="black", 
                                          linetype = "dashed"),
        panel.background = element_rect(fill = alpha("orange", 0.1),
                                        colour = alpha("orange", 0.1),
                                        size = 0.5, linetype = "solid"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
grid.arrange(upper, lower_1, lower_2, widths = c(0.33,1.33),
             layout_matrix = rbind(c(1, 2),
                                   c(1, 3)))


# 5.3 heatmap
annot <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
c3_color <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("geo", "type"))
for (sample in names(filtered[,-1])) {
  type <- annot$cit.coloncancermolecularsubtype[which(annot$geo_accession == sample)]
  if (type == "C3") {
    curr_col <- "red"
  } else {
    curr_col <- "blue"
  }
  c3_color <- rbind(c3_color, data.frame(geo = sample, type = curr_col))
}

heatmap.2(data.matrix(filtered)[,-1], 
          ColSideColors = c3_color$type, 
          trace = "none", 
          margins = c(6,2),
          col = rev(heat.colors(12)),
          cexCol = 0.75,
          xlab = "Sample name", ylab = "Probe set",
          labRow = "",
          key.xlab = "Gene expression level",
          lwid=c(1, 6),
          lhei=c(1, 4))


# 5.4 statistics
# apply this to every row
superT <- function(rows) {
  result <- t.test(rows[,as.vector(which(groups == 1))], 
                   rows[,as.vector(which(groups == 2))])
  return(result)
}

ttest_names <- c("probeID", "t_stat", "p_valu", "p_adju")
ttests <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), ttest_names)
for (rows in c(1:length(filtered$X))) {
  results <- superT(filtered[rows,-1])
  tempDF <- data.frame(probeID = filtered[rows, 1], t_stat = results$statistic,
                       p_valu = results$p.value, 
                       p_adju = p.adjust(results$p.value, 
                                         "fdr", dim(filtered)[2] - 1))
  ttests <- rbind(ttests, tempDF)
}

print(paste("There are", sum(ttests$p_adju < 0.05), "probe sets with p <0.05"))
write.csv(ttests, "clustered_ttests.csv")
