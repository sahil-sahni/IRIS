library(tidyverse)
library(magrittr)
library(plyr)
library(readxl)
library(ModelMetrics)
library(aod)
library(cutpointr)
library(ggpubr)
library(survminer)
###

# Updated from original submission 05/22/24 SS

## ------- Functions -------

survival_tide_quantile <- function(test, name, features){
  score_df = merge(features, test, by.x =1, by.y=1) %>% dplyr::select(1,3,5,6) 
  score_df = score_df %>% mutate(pred=(TIDE>quantile(TIDE,1/2,na.rm=T))*1) %>% # when score > median 
    mutate(time=as.numeric(test[,4]), status=as.numeric(test[,3]))
  return(score_df)
}
calcsurv_low <- function(scoredf, title=NA){
  scoredf = scoredf %>% mutate(strata = ifelse(pred==1,'high','low'))
  fit <- survfit(Surv(time, status) ~ strata, data = scoredf)
  ggtheme = theme_pubr() + theme(plot.title = element_text(hjust = 0.5))
  g=ggsurvplot(fit,pval=TRUE, censor=TRUE, ylab="Survival", legend=c(0.8,0.8), data=scoredf, palette=c('#88BECD','#ffb400'), ggtheme=ggtheme, legend.title='risk group', title=title)+xlab('days')+ylab('survival probability')
  return(g)
}
## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## ------- TIDE OVERALL SURVIVAL -------
## ------- Input -------

## load ICB data sets
test_input = ICB_OS_input[[1]]

## load TIDE
TIDE_score = TIDE_score[-7:-8]

## ------- Calculate Survival -------
surv_output=lapply(1:length(test_input), function(or) survival_tide_quantile(test_input[[or]], names(test_input[or]), TIDE_score[[or]]))
names(surv_output) = names(test_input)

## ------- Supp 3c -------
score_df = do.call(rbind, surv_output[c(2,3,5,6)])
sup3c = calcsurv_low(score_df, 'TIDE OS')
sup3c = plot(sup3c$plot)

## ------- TIDE PFS SURVIVAL -------

## load external data sets
test_input = ICB_PFS_input[[1]]

## load TIDE
TIDE_score = TIDE_score[-7:-8]

## ------- Calculate Survival -------
surv_output=lapply(1:length(test_input), function(or) survival_tide_quantile(test_input[[or]], names(test_input[or]), TIDE_score[[or]]))
names(surv_output) = names(test_input)

## ------- Supp 3d -------
score_df = do.call(rbind, surv_output[c(2,3,5)])  #remove puch (6)
sup3d = calcsurv_low(score_df, 'TIDE PFS')
sup3d = plot(sup3d$plot)

