

ttest_results<-read.csv("/projectnb/bf528/users/hedgehog/project_1/noise_filtering/clustered_ttests.csv")
head(ttest_results)
ttest_results$lfch<-logratio2foldchange(ttest_results$p_adju, base=2)
ttest_r_pvalue<-arrange(ttest_results, desc(ttest_results$p_adju))


