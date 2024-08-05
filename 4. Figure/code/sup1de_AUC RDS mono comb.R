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

## load ICB data sets (mono or comb therapy specific)
test_input = ICB_input[[2]][-13] #remove Auslander mono post (no responder)

## IMPRES results
IMPRES_AUC_supp
IMPRES_AUC_supp[12,1] = 'PUCH_mono_pre'

## TIDE results
TIDE_AUC_supp
TIDE_PUCH =  TIDE_AUC[6,]
TIDE_PUCH[1,1] = 'PUCH_mono_pre'
TIDE_AUC_supp_2 = rbind(TIDE_PUCH,TIDE_AUC_supp)

## ------- Calculate AUC for RDI and RAI (mono or comb therapy specific) -------

RDS_AUC=lapply(1:length(test_input), function(or) calculate_auc(test_input[[or]], names(test_input[or]), feature_RDI_mc[[or]], '<')) #< RDI
RDS_AUC=do.call(rbind, RDS_AUC) 
RDS_AUC[14,1] = 'PUCH_mono_pre'
colnames(RDS_AUC)[2] = 'RDS'

## ------- Combine AUC -------
bar_input = merge(RDS_AUC, IMPRES_AUC_supp, by.x=1, by.y=1, all=T) %>% merge(., TIDE_AUC_supp_2, by.x=1, by.y=1, all=T) %>% reshape2::melt(.) %>% set_colnames(c('cohort', 'method', 'AUC')) %>% 
  mutate(cohort=gsub('_',' ', cohort)) %>% mutate(cohort=gsub(" on", " post", cohort)) %>%
  mutate(cohort=factor(cohort, levels=c('Auslander comb post', 'Gide comb post', 'Riaz comb post',
                                        'Gide mono post', 'Riaz mono post',
                                        'Auslander comb pre', 'Gide comb pre', 'Riaz comb pre', 'Liu comb pre',
                                        'Auslander mono pre', 'Gide mono pre', 'Riaz mono pre', 'Liu mono pre', 'PUCH mono pre'))) 

## ------- Supp 1D: Monotherapy -------
bar_input_mono = bar_input %>% subset(cohort %in% c('Gide mono post', 'Riaz mono post',
                                                    'Auslander mono pre', 'Gide mono pre', 'Riaz mono pre', 'Liu mono pre', 'PUCH mono pre')) 

sup1d=ggplot(data=bar_input_mono, aes(x=cohort, y=AUC, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(), color='black', size=0.3)+
  geom_text(aes(label=round(AUC,2)), vjust=-0.3, size=5/.pt, position = position_dodge(0.9))+ 
  theme_pubr() +
  theme(legend.position = "right", axis.title.x = element_blank(), axis.text.x=element_text(angle=-30,hjust=0))+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.05))+
  scale_fill_manual(values=c('#D0759F','#b3bfd1','#d7e1ee'))
plot(sup1d)

bar_input_mono %>% group_by(method) %>% dplyr::summarise(mean=mean(AUC))
bar_input_mono %>% subset(cohort %in% c('Auslander mono pre', 'Gide mono pre', 'Riaz mono pre', 'Liu mono pre', 'PUCH mono pre')) %>% group_by(method) %>% dplyr::summarise(mean=mean(AUC))

## ------- Supp 1D: Combtherapy -------
bar_input_comb = bar_input %>% subset(cohort %in% c('Auslander comb post', 'Gide comb post', 'Riaz comb post',
                                                    'Auslander comb pre', 'Gide comb pre', 'Riaz comb pre', 'Liu comb pre'))

sup1e=ggplot(data=bar_input_comb, aes(x=cohort, y=AUC, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(), color='black',size=0.3)+
  geom_text(aes(label=round(AUC,2)), vjust=-0.3, size=5/.pt, position = position_dodge(0.9))+ 
  theme_pubr() +
  theme(legend.position = "right", axis.title.x = element_blank(), axis.text.x=element_text(angle=-30,hjust=0))+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.1))+
  scale_fill_manual(values=c('#D0759F','#b3bfd1','#d7e1ee'))
plot(sup1e)

bar_input_comb %>% group_by(method) %>% dplyr::summarise(mean=mean(AUC))
bar_input_comb %>% subset(cohort %in% c('Auslander comb pre', 'Gide comb pre', 'Riaz comb pre', 'Liu comb pre')) %>% group_by(method) %>% dplyr::summarise(mean=mean(AUC))
