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

## load external data sets
test_input = ICB_input[[1]]

## ------- Calculate AUC for RDI -------

## calculate AUC for RDI features
AUC_RDS=lapply(1:length(test_input), function(or) calculate_auc(test_input[[or]], names(test_input[or]), feature_RDI[[or]], '<')) #< RDI
AUC_RDS=do.call(rbind, AUC_RDS) 
AUC_RDS

## ------- Combine AUC -------
LIRICS_output = AUC_RDS
colnames(LIRICS_output)[2] = 'RDS'
bar_input =  merge(LIRICS_output, TIDE_AUC, by.x=1, by.y=1, all=T) %>% merge(., sig_AUC, by.x=1, by.y=1, all=T) %>% reshape2::melt(.) %>% set_colnames(c('cohort', 'method', 'AUC')) %>% 
  mutate(cohort=gsub('_',' ', cohort))%>% mutate(cohort=ifelse(cohort=='Puch pre','PUCH pre', cohort)) %>%  
  mutate(cohort=gsub("on", "post", cohort)) %>%
  mutate(cohort=factor(cohort, levels=c('Auslander post', 'Gide post', 'Riaz post','Auslander pre', 'Gide pre', 'Riaz pre', 'PUCH pre', 'Liu pre'))) %>% 
  subset(method %in% c('RDS','IMPRES','TIDE','MPS','Cytotoxic','resF'))%>%
  mutate(method=factor(method, levels=c('RDS','IMPRES','TIDE','MPS','Cytotoxic','resF'))) 
  
## ------- Figure 2C -------
fig2c=ggplot(data=bar_input, aes(x=cohort, y=AUC, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(), color='black')+
  geom_text(aes(label=round(AUC,2)), vjust=0.5, hjust=-0.4, size=3.5, angle=90, position = position_dodge(0.9))+ 
  scale_y_continuous(expand=c(0,0),limits=c(0,1.10))+
  theme_pubr() +
  theme(legend.position = "right", axis.title.x = element_blank(), axis.text.x=element_text(angle=-30, hjust=0))+
  scale_fill_manual(values=c('#D0759F','#C284A8','#B392B1', '#A5A1BB', '#96AFC4', '#88BECD'))
plot(fig2c)

## paper statistics
bar_input %>% group_by(method) %>% dplyr::summarise(mean=mean(AUC))

