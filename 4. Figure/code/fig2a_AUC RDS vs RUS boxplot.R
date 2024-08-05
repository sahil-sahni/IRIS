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

# Updated from original submission 05/22/24 SS

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

## ------- Calculate AUC for RDI and RUI-------

## calculate AUC for RDI features
AUC_RDS=lapply(1:length(test_input), function(or) calculate_auc(test_input[[or]], names(test_input[or]), feature_RDI[[or]], '<')) #< RDI
AUC_RDS=do.call(rbind, AUC_RDS)

## calculate AUC for RUI features
AUC_RUS=lapply(1:length(test_input), function(or) calculate_auc(test_input[[or]], names(test_input[or]), feature_RUI[[or]], '>')) #> RUI
AUC_RUS=do.call(rbind, AUC_RUS)

## ------- Figure 2A -------
box_input = merge(AUC_RDS, AUC_RUS, by.x=1, by.y=1, all=T) %>% set_colnames(c('cohort', 'RDS', 'RUS')) %>% reshape2::melt(.) %>% set_colnames(c('cohort', 'method', 'AUC'))

fig2a=ggboxplot(box_input, x ="method", y = "AUC",  fill="method", add='jitter') + 
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ylim(0,1)+
  stat_compare_means(ref.group='RDS', paired=T, method='wilcox.test', method.args = list(alternative = "greater"), bracket.size = 0.4, label='p.signif', comparisons = list(c('RDS','RUS')), tip.length = 0.03, size=5, label.y=0.85)+
  scale_fill_manual(values=c('#D0759F','#88BECD'))
plot(fig2a)

## paper statistics
box_input %>% group_by(method) %>% dplyr::summarise(mean=mean(AUC))

#saveRDS(file='fig2a.rds',fig2a)