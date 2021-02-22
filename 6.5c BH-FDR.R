

install.packages(sgof)
library(sgof)
ttest_wBH<-read.csv("/projectnb/bf528/users/hedgehog/project_1/noise_filtering/clustered_ttests.csv")
head(ttest_results)



ttest_wBH$FDR<-p.adjust(ttest_results$p_valu, method = "hochberg")
head(ttest_wBH)
