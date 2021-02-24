library(dplyr)
library(rlang)
library(purrr)
library(stringr)
library(tidyverse)
library(tidyselect)
library(base)

ttest_results<-read.csv("/projectnb/bf528/users/hedgehog/project_1/noise_filtering/clustered_ttests.csv")
head(ttest_results)

str(ttest_results)

ttest_results_B<-ttest_results %>% top_n(10, desc(p_adju))%>%arrange(p_adju)
str(ttest_results_B)
ttest_Bs<-format(ttest_results_B, digits = 4, scientific = 2)

knitr::kable(
  ttest_Bs[],    caption = "Gene Sets by Lowest p-value (adjusted)."
)


ttest_results<-read.csv("/projectnb/bf528/users/hedgehog/project_1/noise_filtering/clustered_ttests.csv")
head(ttest_results)

str(ttest_results)

ttest_results_top<-ttest_results %>% arrange(p_adju)
str(ttest_results_top)
ttest_T<-format(ttest_results_top, digits = 3, nsmall= 20, exp(30), scientific = TRUE)
top_tp<-tail(ttest_T,10)


knitr::kable(
  top_tp[],    
  caption = "Gene Sets by Highest p-value (adjusted)."
)


                        
