library(dplyr)
library(rlang)
library(purrr)
library(stringr)
library(tidyverse)
library(tidyselect)
library(base)

install.packages(sgof)
library(sgof)
ttest_wBH<-read.csv("/projectnb/bf528/users/hedgehog/project_1/noise_filtering/clustered_ttests.csv")
head(ttest_wBH)

ttest_wBHm<-format(ttest_wBH, digits = 3, nsmall= 20, exp(30), scientific = TRUE)
ttest_wBHm$FDR<-p.adjust(ttest_wBHm$p_valu, method = "BH")
FDRn10<-head(ttest_wBHm,10)%>%arrange(FDR)
FDRn10
FDRn3<-head(ttest_wBHm,3)%>%arrange(FDR)
FDRn3
summary(FDRn3)


knitr::kable(
  FDRn3[],    
  caption = "Top 3 Gene Sets (with BH-FDR) Analysis."
)
