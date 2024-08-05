library(tidyverse)
library(magrittr)
library(plyr)
library(readxl)
library(ModelMetrics)
library(aod)
library(glmnet)
library(cutpointr)
library(survival)
library(ggpubr)
library(purrr)
library(pROC)

# Updated from original submission 06/14/24 SS

## ------- Functions -------

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## load TCGA data sets
tcga_input = TCGA_input

## load selected features
AR_LIRICS = feature_RDI
feature_AR = AR_LIRICS[c(2,3,5,6,8)] %>% flatten(.) 

## ------- Calculate RDS score for TCGA SKCM  -------
cci_score = calculate_score_scaled_tcga(tcga_input, 'tcga', feature_AR) 

## ------- Figure 4F -------
box_input = cci_score %>% dplyr::select(temperature_binary,scale) %>% na.omit(.) %>%
  mutate(temperature_binary=mapvalues(temperature_binary, from=c('hot','cold'), to=c('hot niche', 'cold niche'))) %>%
  mutate(temperature_binary = temperature_binary %>% factor(., levels=c('hot niche', 'cold niche')))
fig4f=ggplot(box_input, aes(x=temperature_binary, y=scale, fill=temperature_binary, color=temperature_binary))+
  geom_violin(trim=FALSE, fill='white', show.legend = F)+
  geom_boxplot(width=0.3, color='black', show.legend = F)+
  theme_pubr()+
  xlab(NULL)+
  ylab('RDS')+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  scale_fill_manual(values=c('#D0759F','#88BECD'))+
  scale_color_manual(values=c('#D0759F','#88BECD'))+
  stat_compare_means(ref.group='hot niche', paired=F, method='wilcox.test', method.args = list(alternative = "greater"), bracket.size = 0.4, label='p', 
                     comparisons = list(c('hot niche','cold niche')), tip.length = 0.03, size=6, label.y=3.5)
plot(fig4f)

## ------- Calculate AUC classifying hot vs. cold  -------
hc_score= cci_score %>% mutate(resp_binary=ifelse(temperature_binary=="hot", 1,0))
auc=pROC::auc(as.numeric(hc_score$resp_binary), as.numeric(hc_score$score), direction = '<')

rocobj <- pROC::roc(hc_score$resp_binary, hc_score$score, direction = '<')

## ------- Supp 6F -------
sup6f=ggroc(rocobj,color='#D0759F',size=0.75)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color = "grey50", linetype = "dashed")+ 
  annotate("text", x = 0.25, y = 0.05, label = paste0("AUC = ", round(auc, 2)))+theme_pubr()
plot(sup6f)
