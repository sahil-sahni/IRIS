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
  res=suppressWarnings(suppressMessages(as.numeric(pROC::auc(as.numeric(score_df$response), as.numeric(score_df$score), direction = direction)))) # < RDI; > RUI
  return(data.frame(cohort=name, LIRICS=res))
}

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## load external data sets
test_input = ICB_input[[1]]

## ------- Calculate AUC for RDI and RUI -------

## calculate AUC for RDI features
RDS_AUC=lapply(1:length(test_input), function(or) calculate_auc(test_input[[or]], names(test_input[or]), feature_RDI[[or]], '<')) #< AR
RDS_AUC=do.call(rbind, RDS_AUC)

## calculate AUC for RUI features
RUS_AUC=lapply(1:length(test_input), function(or) calculate_auc(test_input[[or]], names(test_input[or]), feature_RUI[[or]], '>')) #> R
RUS_AUC=do.call(rbind, RUS_AUC)

## ------- Combine AUC -------
bar_input = merge(RDS_AUC, RUS_AUC, by.x=1, by.y=1, all=T) %>% set_colnames(c('cohort', 'RDS', 'RUS')) %>% 
  reshape2::melt(.) %>% set_colnames(c('cohort', 'method', 'AUC')) %>% mutate(cohort=gsub('_',' ', cohort)) %>% 
  mutate(cohort=ifelse(cohort == 'Puch pre','PUCH pre',cohort)) %>%
  mutate(cohort=gsub('_',' ', cohort)) %>% mutate(cohort=gsub(" on", " post", cohort)) %>%
  mutate(cohort=factor(cohort, levels=c('Auslander post', 'Gide post', 'Riaz post',
                                        'Auslander pre', 'Gide pre', 'Riaz pre', 'Liu pre', 'PUCH pre'))) 

## ------- Supp 1A -------
sup1a=ggplot(data=bar_input, aes(x=cohort, y=AUC, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(), color='black')+
  geom_text(aes(label=round(AUC,2)), vjust=-0.3, size=3.5, position = position_dodge(0.9))+ 
  scale_y_continuous(expand=c(0,0),limits=c(0,1.05))+
  theme_pubr() +
  theme(legend.position = "right", axis.title.x = element_blank(), axis.text.x=element_text(angle=-30,hjust=0))+
  scale_fill_manual(values=c('#D0759F','#88BECD')) 
plot(sup1a)

bar_input %>% group_by(method) %>% dplyr::summarise(mean=mean(AUC))


