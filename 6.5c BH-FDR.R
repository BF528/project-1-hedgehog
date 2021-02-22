

install.packages(sgof)
library(sgof)
ttest_results<-read.csv("/projectnb/bf528/users/hedgehog/project_1/noise_filtering/clustered_ttests.csv")
head(ttest_results)
p.adjust(ttest_results$p_valu)




