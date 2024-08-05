library(tidyverse)
library(magrittr)
library(plyr)
library(readxl)
library(ModelMetrics)
library(aod)
library(cutpointr)
library(survival)
library(ggpubr)
library(purrr)
library(pROC)

# Updated from original submission 05/22/24 SS

## ------- Functions -------

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## load TCGA data sets
tcga_input = TCGA_input

## load selected features
AR_LIRICS = feature_RDI
feature_AR = AR_LIRICS[c(2,3,5,6,8)] %>% flatten(.) 

cci_score = calculate_score_scaled_tcga(tcga_input, 'tcga', feature_AR) 

## ------- Supp 6e -------
box_input = cci_score %>% dplyr::select(temperature,scale) %>% na.omit(.) %>%
  mutate(temperature=tolower(temperature)) %>%
  mutate(temperature = temperature %>% factor(., levels=c('non-brisk multifocal', 'non-brisk focal', 'brisk band-like', 'brisk diffuse')))
sup6e=ggplot(box_input, aes(x=temperature, y=scale, fill=temperature, color=temperature))+
  geom_violin(trim=FALSE, fill='white', show.legend = F)+
  geom_boxplot(width=0.20, color='black', show.legend = F)+
  theme_pubr()+
  xlab(NULL)+
  ylab('RDS')+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  theme(axis.text.x = element_text(angle=30,hjust=1))+
  scale_fill_manual(values=c('#88BECD','#75D0C1','#ffb400','#D0759F'))+
  scale_color_manual(values=c('#88BECD','#75D0C1','#ffb400','#D0759F'))+
  stat_compare_means(paired=F, method='wilcox.test', method.args = list(alternative = "greater"), bracket.size = 0.4, label='p.signif', 
                     comparisons = list(c('brisk band-like','non-brisk focal'), 
                                        c('brisk band-like', 'non-brisk multifocal'),
                                        c('brisk diffuse', 'non-brisk focal'),
                                        c('brisk diffuse', 'non-brisk multifocal')), tip.length = 0.03, size=6, label.y=c(4,5.5,7,8.5))
plot(sup6e)


