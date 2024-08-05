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

survival_signature_quantile <- function(test, name, features){
  features=features[[name]]
  score_df = merge(features, test, by.x =1, by.y=1) %>% dplyr::select(1:14)
  for(z in 2:11){
    score_df[,z] = (score_df[,z]>quantile(score_df[,z],1/2,na.rm=T))*1
  }
  score_df = score_df %>%
    mutate(time=as.numeric(.[,14]), status=as.numeric(.[,13]))
  return(score_df)
}

calcsurv_low <- function(scoredf, title=NA){
  scoredf = scoredf %>% mutate(strata = ifelse(pred==1,'high','low'))
  fit <- survfit(Surv(time, status) ~ strata, data = scoredf)
  ggtheme = theme_pubr() + theme(plot.title = element_text(hjust = 0.5))
  g=ggsurvplot(fit,pval=TRUE, censor=TRUE, ylab="Survival", legend=c(0.8,0.8), data=scoredf, palette=c('#88BECD','#ffb400'), ggtheme=ggtheme, legend.title='risk group', title=title)+xlab('days')+ylab('survival probability')
  return(g)
}

calcsurv_high <- function(scoredf, title=NA){
  scoredf = scoredf %>% mutate(strata = ifelse(pred==1,'low','high'))
  fit <- survfit(Surv(time, status) ~ strata, data = scoredf)
  ggtheme = theme_pubr() + theme(plot.title = element_text(hjust = 0.5))
  g=ggsurvplot(fit,pval=TRUE, censor=TRUE, ylab="Survival", legend=c(0.8,0.8), data=scoredf, palette=c('#88BECD','#ffb400'), ggtheme=ggtheme, legend.title='risk group', title=title)+xlab('days')+ylab('survival probability')
  return(g)
}

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## ------- ICB Overall Survival -------
## ------- Input -------

## load ICB OS data sets
test_input = ICB_OS_input[[1]]

## load signatures
signatures = sig_score

## ------- Calculate Score -------
surv_output=lapply(1:length(test_input), function(or) survival_signature_quantile(test_input[[or]], names(test_input[or]), signatures))
names(surv_output) = names(test_input)

## ------- Supp 3EGIKMOQSUW -------
score_df = do.call(rbind, surv_output[c(2,3,5,6)])

fig.letter = c('e','g','i','k','m','o','q','s','u','w')
l=1
for(i in 2:11){
  sel_sig_df =  score_df %>% dplyr::select(1,pred=i,status,time)
  sup = calcsurv_high(sel_sig_df, paste(colnames(score_df)[[i]], 'OS', sep=' '))
  if (colnames(score_df)[[i]] %in% c('MPS','resF')){
    sup = calcsurv_low(sel_sig_df, paste(colnames(score_df)[[i]], 'OS', sep=' '))
  }
  sup = plot(sup$plot)
  filename = paste('sup3',fig.letter[[l]],'.rds',sep='')
  #saveRDS(file=filename, sup)
  l = l+1
}

## ------- ICB PFS Survival -------
## ------- Input -------

## load ICB PFS data sets
test_input = ICB_PFS_input[[1]]

## load signatures
signatures = sig_score

## ------- Calculate Survival -------
surv_output=lapply(1:length(test_input), function(or) survival_signature_quantile(test_input[[or]], names(test_input[or]), signatures))
names(surv_output) = names(test_input)

## ------- Supp 3FHJLNPRTVX -------
score_df = do.call(rbind, surv_output[c(2,3,5)])

fig.letter = c('f','h','j','l','n','p','r','t','v','x')
l=1
for(i in 2:11){
  sel_sig_df =  score_df %>% dplyr::select(1,pred=i,status,time)
  sup = calcsurv_high(sel_sig_df, paste(colnames(score_df)[[i]], 'PFS', sep=' '))
  if (colnames(score_df)[[i]] %in% c('MPS','resF')){
    sup = calcsurv_low(sel_sig_df, paste(colnames(score_df)[[i]], 'PFS', sep=' '))
  }
  sup = plot(sup$plot)
  filename = paste('sup3',fig.letter[[l]],'.rds',sep='')
  #saveRDS(file=filename, sup)
  l = l+1
}
