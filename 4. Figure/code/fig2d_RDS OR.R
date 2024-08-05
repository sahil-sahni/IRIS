library(tidyverse)
library(magrittr)
library(plyr)
library(readxl)
library(ModelMetrics)
library(aod)
library(cutpointr)
library(survival)
library(ggpubr)
library(cutpointr)
###

# Updated from original submission 05/23/24 SS

## ------- Functions -------
calculate_cutoff <- function(train, features){
  scr=lapply(1:length(features), function(x){
    fx = features[[x]]
    tx = train[which(names(train) == names(features[x]))][[1]]
    
    if (length(fx) == 0){
      return(NA)
    }
    
    opt_df = tx[,-1:-2] %>% .[colnames(.) %in% fx] %>% apply(1,function(x){sum(x,na.rm=T)}) %>% data.frame(tx[,1:2], score=./length(fx)) %>% mutate(score=scale(score), response=as.numeric(response))
    opt_scr = cutpointr(opt_df, score, response, direction = ">=", pos_class=1, neg_class=0, method=maximize_metric, metric=sum_sens_spec, na.rm=T)
    return(opt_scr$optimal_cutpoint)
  })
  names(scr) = names(features)
  return(scr)
}

or_calc <- function(test, name, features, cutoff){
  features_name = names(features) %>% unlist(.)
  cutoff_name = names(cutoff) %>% unlist(.)
  
  features_1 = features[which(name == features_name)] %>% flatten(.)
  cutoff_1 = cutoff[which(name == cutoff_name)] %>% flatten(.)
  
  cutoff_m = unlist(cutoff_1) %>% mean(., na.rm=T)

  rank_df=lapply(1:length(features_1), function(r){
    rank = test[,-1:-2] %>% .[colnames(.) %in% features_1[[r]]] %>% apply(1, sum) 
    rank = rank/length(features_1[[r]]) 
    rank=scale(rank)
    return(rank)
  })
  rank_df=do.call(cbind, rank_df)
  score_df = rank_df %>% apply(1,function(x){sum(x,na.rm=T)}) %>% data.frame(response=as.numeric(test[,2]), score=.) %>% mutate(pred=ifelse(score>cutoff_m,1,0)) #/length(features

  tp=nrow(subset(score_df, response == 1 & pred ==1))
  tn=nrow(subset(score_df, response == 0 & pred ==0))
  fp=nrow(subset(score_df, response == 0 & pred ==1))
  fn=nrow(subset(score_df, response == 1 & pred ==0))
  
  res = (tp+tn)/(fn+fp)

  return(data.frame(cohort=name, or=res))
}

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## load ICB data sets
test_input = ICB_input[[1]]
train_input = ICB_input[[1]]


## ------- Calculate odds ratio cut-off -------
cutoff_RDI=lapply(1:length(train_input), function(or) calculate_cutoff(train_input, feature_RDI[[or]]))
names(cutoff_RDI) = names(train_input)

## ------- Calculate odds ratio -------
LIRICS_output=lapply(1:length(test_input), function(or) or_calc(test_input[[or]], names(test_input[or]), feature_RDI, cutoff_RDI))
LIRICS_output=do.call(rbind, LIRICS_output)

## paper statistics
mean(LIRICS_output$or)

## ------- Figure 2D -------
bar_input = LIRICS_output %>% mutate(cohort=gsub('_',' ', cohort)) %>% mutate(cohort=ifelse(cohort=='Puch pre','PUCH pre', cohort)) %>%
  mutate(cohort=gsub("on", "post", cohort)) %>%
  mutate(cohort=factor(cohort, levels=c('Auslander post', 'Gide post', 'Riaz post','Auslander pre', 'Gide pre', 'Riaz pre', 'PUCH pre', 'Liu pre')))
fig2d=ggplot(data=bar_input, aes(x=cohort, y=or, fill=T)) +
  theme_pubr()+
  geom_bar(stat="identity", position=position_dodge(), color='black')+
  geom_text(aes(label=round(or,2)), vjust=-0.3, size=3.5, position = position_dodge(0.9))+ 
  ylab('odds ratio')+  
  scale_y_continuous(expand=c(0,0),limits=c(0,5.5))+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle=-30,hjust=0))+
  scale_fill_manual(values=c('#D0759F')) + 
  geom_hline(yintercept=1, linetype="dashed", color = "black", linewidth=0.5)
plot(fig2d)


