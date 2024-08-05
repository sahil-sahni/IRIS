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
library(ggrepel)
###

# Updated from original submission 05/23/24 SS

## ------- Functions -------
fisher_annotations <- function(selected_features){
  all_annotations = selected_features %>% unlist(.) %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
    set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(GenePair=paste(L,R,sep='_')) %>% merge(., LR_anno, by.x='GenePair', by.y='GenePair')  %>% mutate(annotations=ifelse(is.na(annotations),'None',annotations))
  all_annotations = all_annotations$annotations %>% unique(.) 
  
  ct_auc = lapply(1:length(all_annotations), function(x){
    an = all_annotations[[x]]
    
    an_enriched = selected_features %>% unlist(.) %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
      set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(GenePair=paste(L,R,sep='_')) %>% merge(., LR_anno, by.x='GenePair', by.y='GenePair') %>% mutate(annotations=ifelse(is.na(annotations),'None',annotations)) %>%
      subset(annotations == an) %>% .$GenePair %>% unique(.)
    
    an_lirics = lirics_relevant %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
      set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(GenePair=paste(L,R,sep='_')) %>% merge(., LR_anno, by.x='GenePair', by.y='GenePair') %>% mutate(annotations=ifelse(is.na(annotations),'None',annotations)) %>%
      subset(annotations == an) %>% .$GenePair %>% unique(.) %>% subset(!(. %in% an_enriched))
    
    nan_enriched = selected_features %>% unlist(.) %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
      set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(GenePair=paste(L,R,sep='_')) %>% merge(., LR_anno, by.x='GenePair', by.y='GenePair') %>% mutate(annotations=ifelse(is.na(annotations),'None',annotations)) %>%
      subset(annotations != an) %>% .$GenePair %>% unique(.)
    
    nan_lirics = lirics_relevant %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
      set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(GenePair=paste(L,R,sep='_')) %>% merge(., LR_anno, by.x='GenePair', by.y='GenePair') %>% mutate(annotations=ifelse(is.na(annotations),'None',annotations)) %>%
      subset(annotations != an) %>% .$GenePair %>% unique(.) %>% subset(!(. %in% nan_enriched))
    
    mat = data.frame(
      "enriched" = c(length(an_enriched), length(an_lirics)), 
      "not_enriched" = c(length(nan_enriched), length(nan_lirics)),
      row.names= c("AR", "LIRICS"),
      stringsAsFactors = F
    ) %>% as.matrix(.)
    
    # 4 cell matrix: 
    #         enriched  not enriched
    # RDI      S is 1    L is 0
    # LIRICS  L is 1    S is 0
    res = data.frame(p=fisher.test(mat, alternative='g')$p.value)
    row.names(res)[1] = an
    return(res)
    
    
    
  })
  
  return(ct_auc %>% do.call(rbind,.))
  
}
## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

# load lirics annotations
lirics_annotations = read_excel(paste(PATH,'NIHMS1770717-supplement-2.xlsx',sep=''),sheet=1,skip=1) # download from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8983586/ 

## load ICB data sets
test_input = ICB_input

# load ligand receptor annotations
LR_anno = lirics_annotations
LR_anno = LR_anno %>% mutate(GenePair = paste(ligand, receptor, sep="_"))

## ------- Calculate Annotation Enrichment -------
an_enrichment=fisher_annotations(feature_RDI) %>% mutate(annotations=row.names(.)) %>% arrange(p) %>% 
  mutate(rank=1:nrow(.)) %>% mutate(P=log2((1/p)))%>%
  mutate(annotations=factor(annotations, levels=c('chemotaxis', 'cell-adhesion/LTEM', 'pro-inflammatory', 'activating/co-stimulatory', 'checkpoint/inhibitory', 'None')))%>%
  subset(annotations != 'None') 

## ------- Figure 3C -------
fig3c=ggplot(data=an_enrichment, aes(y=P, x=annotations,fill=annotations)) +
  geom_bar(stat="identity", position=position_dodge(), color='black',width=0.5)+
  geom_hline(yintercept=log2(1/.05), linetype="dashed", color = "black", linewidth=0.5)+ 
  geom_text(aes(label=round(P,1)), vjust=-0.3, size=3.5, position = position_dodge(0.9))+ 
  theme_pubr() +
  theme(legend.position = 'none', axis.title.x = element_blank())+
  scale_fill_manual(values=c("darkorange","aquamarine2","darkorchid","deeppink","deepskyblue"))+
  ylab('-log2(p)')+
  scale_y_continuous(expand=c(0,0),limits=c(0,12))


plot(fig3c)

