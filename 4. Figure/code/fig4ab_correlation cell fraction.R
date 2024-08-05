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

## ------- Functions -------

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## load TCGA data sets
tcga_input = TCGA_input

## load selected features
AR_LIRICS = feature_RDI
feature_AR = AR_LIRICS[c(2,3,5,6,8)] %>% flatten(.) ##CHECK THIS AGAIN

## ------- Calculate RDS score for TCGA SKCM  -------
cci_score = calculate_score_scaled_tcga(tcga_input, 'tcga', feature_AR) 

## ------- Figure 4  -------
#### FIGURE 4A: CD8+ fraction####
scatter_input = cci_score %>% dplyr::select(scale, cd8_frac) %>% na.omit(.)
fig4a=ggscatter(scatter_input, x = "cd8_frac", y = "scale", color='lightgray',
            add = "reg.line",  # Add regressin line
            add.params = list(color = "#D0759F", fill = "pink"), # Customize reg. line
            conf.int = TRUE)+
  theme_pubr()+
  xlab('CD8+ T cell deconvolved fraction')+
  ylab('RDS')+
  stat_cor(method='pearson', p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0, label.y.npc = 0)
plot(fig4a)
cor.test(scatter_input$scale, scatter_input$cd8_frac, method='pearson')$p.value
#### FIGURE 4B: malignant fraction  ####
scatter_input = cci_score %>% dplyr::select(scale, mal_frac) %>% na.omit(.)
fig4b=ggscatter(scatter_input, x = "mal_frac", y = "scale", color='lightgray',
                add = "reg.line",  # Add regressin line
                add.params = list(color = "#88BECD", fill = "light blue"), # Customize reg. line
                conf.int = TRUE)+
  theme_pubr()+
  xlab('malignant deconvolved fraction')+
  ylab('RDS')+
  stat_cor(method='pearson', p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0, label.y.npc = 0)
plot(fig4b)
cor.test(scatter_input$scale, scatter_input$mal_frac, method='pearson')$p.value

