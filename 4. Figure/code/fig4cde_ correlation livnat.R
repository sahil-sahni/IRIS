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

## ------- Functions -------

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

#LIVNAT resF scores
livnat_tcga = readRDS('TCGA_SKCM.RDS') # download data from https://portals.broadinstitute.org/single_cell/study/melanoma-immunotherapy-resistance 

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
livnat_signatures$purity2 = livnat_tcga$purity
livnat_score = merge(livnat_signatures, livnat_score, by.x='sample', by.y='sample') 
livnat_score=livnat_score %>% mutate(TIL=resF-resF.minus.TIL)

## ------- Figure 4  -------
#### FIGURE 4C: TIL####
scatter_input = livnat_score %>% dplyr::select(scale, TIL) %>% na.omit(.)
fig4c=ggscatter(scatter_input, x = "TIL", y = "scale", color='lightgray', size=1,
            add = "reg.line",  # Add regressin line
            add.params = list(color = "#D0759F", fill = "pink", size=0.75), # Customize reg. line
            conf.int = TRUE)+
  theme_pubr()+
  xlab('T cell infiltration (TIL)')+
  ylab('RDS')+
  stat_cor(method='pearson', p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0, label.y.npc = 0, size=7/.pt)
plot(fig4c)
cor.test(scatter_input$scale, scatter_input$TIL, method=c('pearson'))$p.value

#### FIGURE 4D: exc-TIL####
scatter_input = livnat_score %>% dplyr::select(scale, TIL, exc) %>% mutate(exc.minus.TIL= exc-TIL) %>% na.omit(.)
fig4d=ggscatter(scatter_input, x = "exc.minus.TIL", y = "scale", color='lightgray', size=1,
                add = "reg.line",  # Add regressin line
                add.params = list(color = "#88BECD", fill = "light blue", size=0.75), # Customize reg. line
                conf.int = TRUE)+
  theme_pubr()+
  xlab('T cell exclusion program - TIL')+
  ylab('RDS')+
  stat_cor(method='pearson', p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0, label.y.npc = 0, size=7/.pt)
plot(fig4d)
cor.test(scatter_input$scale, scatter_input$exc.minus.TIL, method=c('pearson'))$p.value

#### FIGURE 4E: res-TIL####
scatter_input = livnat_score %>% dplyr::select(scale, TIL, res) %>% mutate(res.minus.TIL=res-TIL) %>% na.omit(.)
fig4e=ggscatter(scatter_input, x = "res.minus.TIL", y = "scale", color='lightgray', size=1,
                add = "reg.line",  # Add regressin line
                add.params = list(color = "#88BECD", fill = "light blue", size=0.75), # Customize reg. line
                conf.int = TRUE)+
  theme_pubr()+
  xlab('post-resistance program - TIL')+
  ylab('RDS')+
  stat_cor(method='pearson', p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0, label.y.npc = 0, size=7/.pt)
plot(fig4e)
cor.test(scatter_input$scale, scatter_input$res.minus.TIL, method=c('pearson'))$p.value
