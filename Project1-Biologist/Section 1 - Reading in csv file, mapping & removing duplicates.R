


raw_data_1<-read.csv("/projectnb/bf528/users/hedgehog/project_1/biologist role/#1 Differential Expression & hgu133plus2.db/hedgehog1_data.csv")




Synced_2GS<- AnnotationDbi::select(hgu133plus2.db, keys = raw_data_1$X, columns = ("SYMBOL"))

Synced_2GS_nodup<-Synced_2GS[!duplicated(Synced_2GS),]
