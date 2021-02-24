library(cluster)
library(gplots)
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
groups <- cutree(fit, k = 2)  # cut tree into 2 clusters
print(table(groups)) # number of samples in each group
plot(fit, hang = -1, cex = 0.8, xlab = "Sample GEO Accession")
rect.hclust(fit, k = 2, border = "red")

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

# heatmap(data.matrix(filtered)[,-1], ColSideColors = c3_color$type, )
heatmap.2(data.matrix(filtered)[,-1], ColSideColors = c3_color$type, 
          trace = "none", col = rev(heat.colors(12)))

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
