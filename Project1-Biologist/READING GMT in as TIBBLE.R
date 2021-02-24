read_hall<-read.gmt("/usr4/bf528/bcole99/PROJECT 1/h.all.v7.2.symbols.gmt")
read_kegg<-read.gmt("/usr4/bf528/bcole99/PROJECT 1/c2.cp.kegg.v7.2.symbols.gmt")
read_go<-read.gmt("/usr4/bf528/bcole99/PROJECT 1/c5.go.v7.2.symbols.gmt")
all_genes<-rbind(read_go,read_hall,read_kegg)
