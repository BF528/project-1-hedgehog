

ttest_results<-read.csv("/projectnb/bf528/users/hedgehog/project_1/noise_filtering/clustered_ttests.csv")
head(ttest_results)

str(ttest_results)

ttest_results_B<-ttest_results %>% top_n(10, desc(p_adju))%>%arrange(p_adju)
str(ttest_results_B)





                        
