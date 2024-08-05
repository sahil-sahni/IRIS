library(tidyverse)
library(magrittr)
library(plyr)
library(readxl)
library(ModelMetrics)
library(aod)
library(cutpointr)
library(survival)
library(ggpubr)
library(ggsignif)
library(rstatix)
library(scales)
###

# Updated from original submission 05/23/24 SS

## ------- Functions -------
patient_score <- function(test, name, features){
  score_df = calculate_score_scaled(test,name,features)
  score_df = score_df  %>% mutate(cohort=name,response=factor(mapvalues(response,from=c(0,1), to=c('NR','R')), levels=c('R','NR'))) 
  return(score_df)
}

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## load ICB data sets
test_input = ICB_input[[1]]

## ------- Calculate Patient Score -------
cohort_score=lapply(1:length(test_input), function(or) patient_score(test_input[[or]], names(test_input[or]), feature_RDI[[or]]))
cohort_score=do.call(rbind,cohort_score) %>% group_by(cohort)

## ------- Figure 2B -------
box_input = cohort_score  %>% mutate(cohort=gsub('_',' ', cohort)) %>% mutate(cohort=ifelse(cohort=='Puch pre','PUCH pre', cohort)) %>% mutate(cohort=gsub("on", "post", cohort)) %>%  mutate(cohort=factor(cohort, levels=c('Auslander post', 'Gide post', 'Riaz post','Auslander pre', 'Gide pre', 'Riaz pre', 'PUCH pre', 'Liu pre'))) 
stat.test = box_input %>% group_by(cohort) %>% wilcox_test(scale ~ response, ref.group = 'R', alternative = 'greater')
stat.test <- stat.test %>%
  add_xy_position(x = "cohort", dodge = 0.8) %>% add_significance(p.col='p')
fig2b = ggboxplot(box_input, x='cohort', y='scale', fill='response') + 
  theme(legend.position = 'right', axis.title.x = element_blank(), axis.text.x = element_text(angle=-30,hjust=0))
fig2b = fig2b + stat_pvalue_manual(stat.test,  label = "p", tip.length = .02, y.position=3.8, bracket.size = 0.4, size=5)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels=label_number(accuracy=0.01))+
  scale_fill_manual(values=c('#D0759F','#88BECD'))+
  ylab('RDS')
plot(fig2b)
