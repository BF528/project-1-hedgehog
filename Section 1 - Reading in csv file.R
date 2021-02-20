


raw_data_1<-read.csv("/projectnb/bf528/users/hedgehog/project_1/biologist role/#1 Differential Expression & hgu133plus2.db/hedgehog1_data.csv")


head(raw_data_1)
raw_data<-read.csv("/project/bf528/project_1/data/differential_expression_results.csv")


Synced_2GS<- AnnotationDbi::select(hgu133plus2.db, keys = raw_data_1$X, columns = ("SYMBOL"))
