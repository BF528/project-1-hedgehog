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
read_hall<-read.gmt("/usr4/bf528/bcole99/PROJECT 1/h.all.v7.2.symbols.gmt")