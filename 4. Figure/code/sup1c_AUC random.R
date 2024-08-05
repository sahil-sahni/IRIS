library(tidyverse)
library(magrittr)
library(plyr)
library(readxl)
library(ModelMetrics)
library(aod)
library(cutpointr)
library(survival)
library(ggpubr)
library(ggplot2)
library(viridis)
library(ggsignif)
###

# Updated from original submission 05/23/24 SS

## ------- Functions -------

calculate_auc <- function(test, name, features, direction){
  score_df = calculate_score(test, name, features)
  res=suppressWarnings(suppressMessages(as.numeric(pROC::auc(as.numeric(score_df$response), as.numeric(score_df$score), direction = direction)))) # < RDI; > RAI
  
  return(data.frame(cohort=name, LIRICS=res))
}

auc_permute_random <- function(test, name, features, direction){ # calculate AUC when permuting background on LIRICS relevant interactions only
  set.seed(20892) #set seed
  
  input_name = strsplit(name, '_') %>% unlist(.) %>% .[[1]]
  features_name = strsplit(names(features),'_') %>% sapply(.,'[[',1) %>% unlist(.)
  
  train = features[which(input_name != features_name)] #ensure ensembl model is not training and testing with the same cohort

  ensembl_df=lapply(1:length(train), function(r){ #calculate score for each ensembl model
    rank = test[,-1:-2]  %>% .[sample(nrow(.))] %>% .[colnames(.) %in% train[[r]]]  %>% apply(1, sum) 
    rank = rank/length(train[[r]])
    return(rank)
  })
  ensembl_df=do.call(cbind, ensembl_df)
  
  score_df = ensembl_df %>% apply(1,function(x){sum(x,na.rm=T)}) %>% data.frame(response=test[,2], score=.)
  res=suppressWarnings(suppressMessages(as.numeric(pROC::auc(as.numeric(score_df$response), as.numeric(score_df$score), direction = direction)))) # < RDI; > RAI
  
  return(data.frame(cohort=name, random=res))
}

auc_response_random <- function(test, name, features, direction){ # calculate AUC when permuting ICB response only

  score_df = calculate_score(test, name, features)
  score_df = score_df %>% mutate(response=sample(test[,2]))
  res=suppressWarnings(suppressMessages(as.numeric(pROC::auc(as.numeric(score_df$response), as.numeric(score_df$score), direction = direction)))) # < RDI; > RAI
  
  return(data.frame(cohort=name, LIRICS=res))
}

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## load ICB data sets
test_input = ICB_input[[1]]

## ------- Calculate AUC -------

## calculate AUC for RDI features
RDS_AUC=lapply(1:length(test_input), function(or) calculate_auc(test_input[[or]], names(test_input[or]), feature_RDI[[or]], '<')) #< RDI
RDS_AUC=do.call(rbind, RDS_AUC) 

## calculate random AUC by permutation
LIRICS_auc_permute=lapply(1:length(test_input), function(or) auc_permute_random(test_input[[or]], names(test_input[or]), feature_RDI[[or]], '<')) #< RDI
LIRICS_auc_permute=do.call(rbind, LIRICS_auc_permute) 

## calculate random AUC by shuffling response
LIRICS_auc_response=lapply(1:length(test_input), function(or) auc_response_random(test_input[[or]], names(test_input[or]), feature_RDI[[or]], '<')) #< RDI
LIRICS_auc_response=do.call(rbind, LIRICS_auc_response)

## ------- Combine AUC -------

bar_input = merge(RDS_AUC, LIRICS_auc_permute, by.x=1, by.y=1, all=T)  %>% merge(., LIRICS_auc_response, by.x=1, by.y=1, all=T) %>% 
  set_colnames(c('cohort', 'RDS', 'permute background', 'permute response')) %>% reshape2::melt(.) %>% set_colnames(c('cohort', 'method', 'AUC'))%>% 
  mutate(cohort=gsub('_',' ', cohort))%>% mutate(cohort=ifelse(cohort=='Puch pre','PUCH pre', cohort)) %>%  
  mutate(cohort=gsub('_',' ', cohort)) %>% mutate(cohort=gsub(" on", " post", cohort)) %>%
  mutate(cohort=factor(cohort, levels=c('Auslander post', 'Gide post', 'Riaz post','Auslander pre', 'Gide pre', 'Riaz pre', 'PUCH pre', 'Liu pre'))) 

## ------- Supp 1C -------

sup1c=ggplot(data=bar_input, aes(x=cohort, y=AUC, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(), color='black')+
  geom_text(aes(label=round(AUC,2)), vjust=-0.3, size=3.5, position = position_dodge(0.9))+ 
  theme_pubr() +
  theme(legend.position = "right", axis.title.x = element_blank(), axis.text.x=element_text(angle=-30, hjust=0))+
  scale_fill_manual(values=c('#D0759F','#b3bfd1','#d7e1ee'))+
  scale_y_continuous(expand=c(0,0),limits=c(0,1))
plot(sup1c)

## paper statistics
bar_input %>% group_by(method) %>% dplyr::summarise(mean=mean(AUC))
