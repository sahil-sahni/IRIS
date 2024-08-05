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
  ggtheme = theme_pubr() + theme(plot.title = element_text(hjust = 0.5))
  g=ggsurvplot(fit,pval=TRUE, censor=TRUE, ylab="Survival", legend=c(0.8,0.8), data=scoredf, palette=c('#88BECD','#D0759F'), ggtheme=ggtheme, legend.title='risk group', title=title)+xlab('days')+ylab('survival probability')
  return(g)
}

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## ------- ICB Overall Survival -------
## ------- Input -------

## load external ICB OS data sets
test_input = ICB_OS_input[[1]]


## ------- Calculate Survival -------
surv_output=lapply(1:length(test_input), function(or) survival_test_quantile(test_input[[or]], names(test_input[or]), feature_RDI[[or]]))
names(surv_output) = names(test_input)

## ------- SUPP 4 (OS) -------
#### SUPP 4A: pre overall survival (Gide pre) ####
score_df = surv_output[[2]]
sup4a = calcsurv(score_df, 'Gide pre OS')
sup4a = plot(sup4a$plot)


#### SUPP 4B: pre overall survival (Liu pre) ####
score_df = surv_output[[3]]
sup4b = calcsurv(score_df, 'Liu pre OS')
sup4b = plot(sup4b$plot)


#### SUPP 4C: pre overall survival (Riaz pre) ####
score_df = surv_output[[5]]
sup4c = calcsurv(score_df, 'Riaz pre OS')
sup4c = plot(sup4c$plot)


#### SUPP 4D: pre overall survival (PUCH pre) ####
score_df = surv_output[[6]]
sup4d = calcsurv(score_df, 'PUCH pre OS')
sup4d = plot(sup4d$plot)


#### SUPP 4E: post overall survival (merged) ####
score_df = do.call(rbind, surv_output[c(1,4)])
sup4e = calcsurv(score_df, 'ICB post OS')
sup4e = plot(sup4e$plot)


#### SUPP 4F: pre overall survival (Gide post) ####
score_df = surv_output[[1]]
sup4f = calcsurv(score_df, 'Gide post OS')
sup4f = plot(sup4f$plot)


#### SUPP 4G: pre overall survival (Riaz post) ####
score_df = surv_output[[4]]
sup4g = calcsurv(score_df, 'Riaz post OS')
sup4g = plot(sup4g$plot)


## ------- ICB PFS SURVIVAL -------

## ------- Input -------
test_input = ICB_PFS_input[[1]]

## ------- Calculate Survival -------
surv_output=lapply(1:length(test_input), function(or) survival_test_quantile(test_input[[or]], names(test_input[or]), feature_RDI[[or]]))
names(surv_output) = names(test_input)

## ------- SUPP 4 (PFS) -------
#### SUPP 4H: pre progression free survival (Gide pre) ####
score_df = surv_output[[2]]
sup4h = calcsurv(score_df, 'Gide pre PFS')
sup4h = plot(sup4h$plot)


#### SUPP 4I: pre progression free survival (Liu pre) ####
score_df = surv_output[[3]]
sup4i = calcsurv(score_df, 'Liu pre PFS')
sup4i = plot(sup4i$plot)


#### SUPP 4J: pre progression free survival (Riaz pre) ####
score_df = surv_output[[5]]
sup4j = calcsurv(score_df, 'Riaz pre PFS')
sup4j = plot(sup4j$plot)


#### SUPP 4K: post progression free survival (merged) ####
score_df = do.call(rbind, surv_output[c(1,4)])
sup4k = calcsurv(score_df, 'ICB post PFS')
sup4k = plot(sup4k$plot)


#### SUPP 4L: pre progression free survival (Gide post) ####
score_df = surv_output[[1]]
sup4l = calcsurv(score_df, 'Gide post PFS')
sup4l = plot(sup4l$plot)


#### SUPP 4M: pre progression free survival (Riaz post) ####
score_df = surv_output[[4]]
sup4m = calcsurv(score_df, 'Riaz post PFS')
sup4m = plot(sup4m$plot)


