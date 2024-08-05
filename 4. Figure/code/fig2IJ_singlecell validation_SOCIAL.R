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
library(ggplot2)
library(pROC)
library(scales)
library(abind)
###

# Updated from original submission 07/02/24 SS

## ------- Functions -------
patient_score <- function(test, name, features){
  score_df = calculate_score_scaled(test, name, features)
  score_df = score_df %>% mutate(cohort=name,response=factor(mapvalues(test[,2],from=c('NAV','RES'), to=c('naïve','post-ICB resistant')), levels=c('naïve','post-ICB resistant'))) 
  
  return(score_df)
}

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## load external data sets
model_input = ICB_sc_input

## load selected features
AR_LIRICS = feature_RDI
AR_LIRICS = AR_LIRICS[c(2,3,5,6,8)] %>% flatten(.)

## ------- Run SOCIAL step 4 -------
source('~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/3. SOCIAL & SPECIAL/SOCIAL & SPECIAL Function/SOCIAL_SPECIAL.R')

test_input = SOCIAL.infer_activity(model_input)
test_input = data.frame(sample=row.names(test_input), response=substr(row.names(test_input),1,3)) %>% cbind(., test_input)


## ------- Calculate Patient Score -------
sc_score=patient_score(test_input, 'livnat_paper', AR_LIRICS)

## ------- Figure 2I -------
box_input = sc_score
fig2i=ggplot(box_input, aes(x=response, y=scale, fill=response, color=response))+
  geom_violin(trim=FALSE, fill=NA, show.legend = F)+
  geom_boxplot(width=0.26, color='black', show.legend = F)+
  theme_pubr()+
  xlab(NULL)+
  ylab('RDS')+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels=label_number(accuracy=0.01))+
  scale_fill_manual(values=c('#D0759F','#88BECD'))+
  scale_color_manual(values=c('#D0759F','#88BECD'))+
  ggtitle('melanoma ICB single cell')+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_compare_means(ref.group='naïve', paired=F, method='wilcox.test', method.args = list(alternative = "greater"), bracket.size = 0.4, label='p', 
                     comparisons = list(c('naïve','post-ICB resistant')), tip.length = 0.03, size=6, label.y=4.5)
plot(fig2i)

#pvalue
box_input %>% wilcox_test(scale ~ response, alternative='greater',ref.group='naïve')

## ------- Calculate AUC classifying Naive vs post-ICB resistant -------
sc_score = sc_score %>% mutate(resp_binary=ifelse(response=="naïve", 1,0))
auc=pROC::auc(as.numeric(sc_score$resp_binary), as.numeric(sc_score$score), direction = '<')

## ------- Figure 2J -------
rocobj <- pROC::roc(sc_score$resp_binary, sc_score$score, direction = '<')
fig2j=ggroc(rocobj,color='#D0759F',size=0.75)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color = "grey50", linetype = "dashed")+ 
  annotate("text", x = 0.25, y = 0.05, label = paste0("AUC = ", round(auc, 2)))+
  theme_pubr()+
  ggtitle('melanoma ICB single cell')+
  theme(plot.title = element_text(hjust = 0.5))
  
plot(fig2j)





