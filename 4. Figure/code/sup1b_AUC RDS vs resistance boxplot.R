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

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## load ICB data sets
test_input = ICB_input[[1]]

## calculate AUC for RDI features
RDS_AUC=lapply(1:length(test_input), function(or) calculate_auc(test_input[[or]], names(test_input[or]), feature_RDI[[or]], '<')) #< RDI
RDS_AUC=do.call(rbind, RDS_AUC) 

## ------- Combine AUC -------
LIRICS_output = RDS_AUC
colnames(LIRICS_output)[2] = 'RDS'
box_input = merge(LIRICS_output, TIDE_AUC, by.x=1, by.y=1, all=T) %>% merge(., sig_AUC, by.x=1, by.y=1, all=T) %>% reshape2::melt(.) %>% set_colnames(c('cohort', 'method', 'AUC')) %>% 
  mutate(cohort=gsub('_',' ', cohort)) %>% 
  mutate(cohort=gsub('_',' ', cohort)) %>% mutate(cohort=gsub(" on", " post", cohort)) %>%
  mutate(cohort=factor(cohort, levels=c('Auslander post', 'Gide post', 'Riaz post','Auslander pre', 'Gide pre', 'Riaz pre', 'Puch pre', 'Liu pre'))) 

## ------- Supp 1B -------
sup1b=ggplot(data=box_input, aes(x=method, y=AUC, fill=method)) +
  geom_boxplot()+geom_jitter()+
  scale_y_continuous(expand=c(0,0),limits=c(-.01,1.1))+
  theme_pubr() +
  theme(legend.position = "right", axis.title.x = element_blank(), axis.text.x=element_text(angle=-30, hjust=0))+
  scale_fill_manual(values=c('#D0759F','#d7e1ee','#d7e1ee','#d7e1ee','#d7e1ee','#d7e1ee','#d7e1ee','#d7e1ee','#d7e1ee','#d7e1ee','#d7e1ee','#d7e1ee'))#+
  stat_compare_means(paired=T, method='wilcox.test', method.args = list(alternative = 'greater'), bracket.size = 0.4, label='p', 
                     comparisons = list(c('RDS','IMPRES'), 
                                        c('RDS', 'TIDE'),
                                        c('RDS', 'PD1'),
                                        c('RDS', 'PDL1'),
                                        c('RDS','CTLA4'),
                                        c('RDS','Tex')), tip.length = 0.03, size=6, label.y=c(1,1.2,1.4,1.6,1.8,2.0))
  

plot(sup1b)

## paper statistics
box_input %>% group_by(method) %>% dplyr::summarise(mean=mean(AUC))

ggplot(data=box_input, aes(x=method, y=AUC, fill=method)) +
  geom_boxplot()+geom_jitter()+
  stat_compare_means(paired=T, method='wilcox.test', method.args = list(alternative = 'greater'), bracket.size = 0.4, label='p', 
                     comparisons = list(c('RDS','TIDE'), 
                                        c('RDS', 'IMPRES'),
                                        c('RDS', 'MPS'),
                                        c('RDS', 'Cytotoxic'),
                                        c('RDS','resF'),
                                        c('RDS','Tin')), tip.length = 0.03, size=6,)
