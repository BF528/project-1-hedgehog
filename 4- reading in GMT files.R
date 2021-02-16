KEGG_1 <- getGmt("/usr4/bf528/bcole99/PROJECT 1/c2.cp.kegg.v7.2.symbols.gmt",
                 collectionType=BroadCollection(category="c3"),
                 geneIdType=SymbolIdentifier())


GO_1 <- getGmt("/usr4/bf528/bcole99/PROJECT 1/c5.go.v7.2.symbols.gmt",
               collectionType=BroadCollection(category="c3"),
               geneIdType=SymbolIdentifier())

HALL_1 <- getGmt("/usr4/bf528/bcole99/PROJECT 1/h.all.v7.2.symbols.gmt",
                 collectionType=BroadCollection(category="c3"),
                 geneIdType=SymbolIdentifier())

head(KEGG_1)
head(GO_1)
head(HALL_1)


