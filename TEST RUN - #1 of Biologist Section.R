raw_data<-read.csv("/project/bf528/project_1/data/differential_expression_results.csv")


Synced_2GS_test <- AnnotationDbi::select(hgu133plus2.db, keys = as.character(row.names(raw_data)), columns = ("SYMBOL"))
head(raw_data)
head(Synced_2GS_test)