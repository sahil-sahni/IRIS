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

## ------- Calculate RDS score for TCGA SKCM  -------
cci_score = calculate_score_scaled_tcga(tcga_input, 'tcga', feature_AR) 

## ------- Load/Calculate Livnat Signature Scores in TCGA  -------
livnat_tme=livnat_tcga$tme %>% as.data.frame(.) %>% mutate(sample=livnat_tcga$samples) %>% mutate(sample=str_replace_all(sample,'\\.','-')) 
livnat_score = merge(cci_score, livnat_tme, by.x='sample', by.y='sample')

livnat_signatures = livnat_tcga$res %>% as.data.frame(.) %>% mutate(sample=livnat_tcga$samples) %>% mutate(sample=str_replace_all(sample,'\\.','-')) 
livnat_score = merge(livnat_signatures, livnat_score, by.x='sample', by.y='sample') 
livnat_score=livnat_score %>% mutate(TIL=resF-resF.minus.TIL)

## ------- Supp 6  -------
#### Supp 6A: resF-TIL ####
scatter_input = livnat_score %>% dplyr::select(scale, resF.minus.TIL) %>% na.omit(.)
sup6a=ggscatter(scatter_input, x = "resF.minus.TIL", y = "scale", color='lightgray',
                add = "reg.line",  # Add regressin line
                add.params = list(color = "#88BECD", fill = "light blue"), # Customize reg. line
                conf.int = TRUE)+
  theme_pubr()+
  xlab('immune resistance program - TIL')+
  ylab('RDS')+
  stat_cor(method='pearson', p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0, label.y.npc = 0)
plot(sup6a)
## paper results
cor.test(scatter_input$scale,scatter_input$resF.minus.TIL,method = 'pearson')$p.value

#### Supp 6B: exc wout TIL ####
scatter_input = livnat_score %>% dplyr::select(scale, exc) %>% na.omit(.)
sup6b=ggscatter(scatter_input, x = "exc", y = "scale", color='lightgray',
            add = "reg.line",  # Add regressin line
            add.params = list(color = "#88BECD", fill = "light blue"), # Customize reg. line
            conf.int = TRUE)+
  theme_pubr()+
  xlab('T cell exclusion program')+
  ylab('RDS')+
  stat_cor(method='pearson', p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0, label.y.npc = 0)
plot(sup6b)
cor.test(scatter_input$scale,scatter_input$exc,method = 'pearson')$p.value

#### Supp 6C: post-res wout TIL ####
scatter_input = livnat_score %>% dplyr::select(scale, res) %>% na.omit(.)
sup6c=ggscatter(scatter_input, x = "res", y = "scale", color='lightgray',
                add = "reg.line",  # Add regressin line
                add.params = list(color = "#88BECD", fill = "light blue"), # Customize reg. line
                conf.int = TRUE)+
  theme_pubr()+
  xlab('post-resistance program')+
  ylab('RDS')+
  stat_cor(method='pearson', p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0, label.y.npc = 0)
plot(sup6c)
cor.test(scatter_input$scale,scatter_input$res,method = 'pearson')$p.value

#### Supp 6D: resF wout TIL ####
scatter_input = livnat_score %>% dplyr::select(scale, resF) %>% na.omit(.)
sup6d=ggscatter(scatter_input, x = "resF", y = "scale", color='lightgray',
                add = "reg.line",  # Add regressin line
                add.params = list(color = "#88BECD", fill = "light blue"), # Customize reg. line
                conf.int = TRUE)+
  theme_pubr()+
  xlab('immune resistance program')+
  ylab('RDS')+
  stat_cor(method='pearson', p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0, label.y.npc = 0)
plot(sup6d)
cor.test(scatter_input$scale,scatter_input$resF,method = 'pearson')$p.value

