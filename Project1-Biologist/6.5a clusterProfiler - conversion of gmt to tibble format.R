#Installs and loads BioManager from BioConductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
library(BiocManager)

#Installs and loads GSEAmining from BioConductor (based on Geneclustering)
BiocManager::install("GSEAmining")
library(GSEAmining)


#Installs and loads clusterProfiler, part of BiocManager in BioConducter. Built in function for converting file to table. (i.e. read.gmt)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
library(clusterProfiler)

read_hall<-read.gmt("/usr4/bf528/bcole99/PROJECT 1/h.all.v7.2.symbols.gmt")
read_go<-read.gmt("/usr4/bf528/bcole99/PROJECT 1/c5.go.v7.2.symbols.gmt")
read_kegg<-read.gmt("/usr4/bf528/bcole99/PROJECT 1/c2.cp.kegg.v7.2.symbols.gmt")
all_genes<-rbind(read_hall,read_go,read_kegg)
dim(all_genes)

head(all_genes)
