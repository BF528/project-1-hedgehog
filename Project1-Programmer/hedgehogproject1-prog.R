#Mackenzie Knox
#2/14/21
#BF528 Project 1

#R version 3.5.1

#install packages
install.packages("BiocManager")
BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::install("sva")
BiocManager::install("AnnotationDbi")
BiocManager::install("hgu133plus2.db")

#load packages
library(tidyverse)
library(ggplot2)
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)


#import CEL files
ds <- ReadAffy(celfile.path = "/projectnb2/bf528/users/hedgehog/project_1/samples_2")
nds <- rma(ds)


dsfit <- fitPLM(ds, normalize=TRUE, background=TRUE)


#Relative Log Expression plot
stats <- RLE(dsfit, type="stats")

df_stats <- as.data.frame(t(stats))
plot_rle <- ggplot(df_stats, aes(x=median)) +
  geom_histogram(color="white", fill="blue")+
  labs(title="Median RLE Scores",
       x="RLE Score", 
       y = "Count")
plot_rle


#Normalized Unscaled Standard Error
nstats <- NUSE(dsfit, type="stats")

df_nstats <- as.data.frame(t(nstats))
plot_nuse <- ggplot(df_nstats, aes(x=median)) +
  geom_histogram(color="white", fill="red")+
  labs(title="Median NUSE Scores",
       x="NUSE Score", 
       y = "Count")
plot_nuse

#import annotation files, correct for batch effects
metadata <- read_csv("/project/bf528/project_1/doc/proj_metadata.csv")
mod1 <- metadata$normalizationcombatmod
mod2 <- model.matrix(~as.factor(mod1))
b <- metadata$normalizationcombatbatch

eds <- exprs(nds)

c <- ComBat(dat=eds, batch=b, mod=mod2)

write.csv(c, file="/projectnb2/bf528/users/hedgehog/project_1/hedgehog1_data.csv")


#Principal Component Analysis
usds <- t(exprs(nds)) #transpose now-matrix
susds<- scale(usds) #scale
sds <- t(susds) #re-transpose

pds <- prcomp(sds, scale. = F, center = F)
pdsrot <- as.data.frame(pds$rotation)


#PC1 v PC2
ggplot(pdsrot, aes(x=PC1, y=PC2))  + 
  geom_point(size=2) +
  labs(title="PC1 v PC2")
summary(pds)

