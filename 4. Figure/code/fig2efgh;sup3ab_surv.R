library(tidyverse)
library(magrittr)
library(plyr)
library(readxl)
library(ModelMetrics)
library(aod)
library(cutpointr)
library(survival)
library(ggpubr)
library(survminer)
###

# Updated from original submission 05/22/24 SS

## ------- Functions -------

survival_test_quantile <- function(test, name, features){
  score_df = calculate_score_surv(test, name, features)
  
  score_df = score_df %>% mutate(pred=(score>quantile(score,1/2,na.rm=T))*1) %>% # when score > median 
    mutate(time=as.numeric(test[,4]), status=as.numeric(test[,3]))
  return(score_df)
  
}

calcsurv <- function(scoredf, title=NA){
  scoredf = scoredf %>% mutate(strata = ifelse(pred==1,'low','high'))
  fit <- survfit(Surv(time, status) ~ strata, data = scoredf)
  print(surv_pvalue(fit, data=scoredf)$pval)
  ggtheme = theme_pubr() + theme(plot.title = element_text(hjust = 0.5))
  g=ggsurvplot(fit,pval=TRUE, censor=TRUE, ylab="Survival", legend=c(0.8,0.8), data=scoredf, palette=c('#88BECD','#D0759F'), ggtheme=ggtheme, legend.title='risk group', title=title)+xlab('days')+ylab('survival probability')
  return(g)
}

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## ------- ################ -------
## ------- ICB OVERALL SURVIVAL (MERGE) -------
## ------- ################ -------

## ------- Input -------

## load external ICB data sets
test_input = ICB_OS_input[[1]]

## ------- Calculate Survival -------
surv_output=lapply(1:length(test_input), function(or) survival_test_quantile(test_input[[or]], names(test_input[or]), feature_RDI[[or]]))
names(surv_output) = names(test_input)

## ------- Figure 2E -------
score_df = do.call(rbind, surv_output[c(2,3,5,6)])
fig2e = calcsurv(score_df, 'pre-treatment ICB OS')
fig2e = plot(fig2e$plot)

## ------- Supp 3a -------
sup3a = calcsurv(score_df, 'RDS OS')
sup3a = plot(sup3a$plot)

## ------- ################ -------
## ------- ICB PROGRESSION FREE SURVIVAL (MERGE) -------
## ------- ################ -------

## ------- Input -------

## load external data sets
test_input = ICB_PFS_input[[1]]
 
## ------- Calculate Survival -------
surv_output=lapply(1:length(test_input), function(or) survival_test_quantile(test_input[[or]], names(test_input[or]), feature_RDI[[or]]))
names(surv_output) = names(test_input)

## ------- Figure 2F -------
score_df = do.call(rbind, surv_output[c(2,3,5)]) #6 is puch remove
fig2f = calcsurv(score_df, 'pre-treatment ICB PFS')
fig2f = plot(fig2f$plot)


## ------- Supp 3b -------
sup3b = calcsurv(score_df, 'RDS PFS')
sup3b = plot(sup3b$plot)


## ------- ################ -------
## ------- TCGA OVERALL SURVIVAL -------
## ------- ################ -------

## ------- Input -------

## load external TCGA data sets
test_input =  TCGA_surv_input[1]

## load selected features
feature_AR = feature_RDI[c(2,3,5,6,8)] %>% flatten(.)

## ------- Calculate Survival -------
surv_output=survival_test_quantile(test_input[[1]], names(test_input[1]), feature_AR)

## ------- Figure 2G -------
fig2g = calcsurv(surv_output, 'TCGA SKCM OS')
fig2g = plot(fig2g$plot)


## ------- ################ -------
## ------- TCGA PROGRESSION FREE SURVIVAL -------
## ------- ################ -------

## ------- Input -------

## load external TCGA data sets
test_input =  TCGA_surv_input[2]

## load selected features
feature_AR = feature_RDI[c(2,3,5,6,8)] %>% flatten(.)

## ------- Calculate Survival -------
surv_output=survival_test_quantile(test_input[[1]], names(test_input[1]), feature_AR)

## ------- Figure 2H -------
fig2h = calcsurv(surv_output, 'TCGA SKCM PFS')
fig2h = plot(fig2h$plot)



