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
fisher_Lcell <- function(selected_features){
  all_cell = selected_features %>% unlist(.) %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
    set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(GenePair=paste(L,R,sep='_')) 
  all_cell = c(all_cell$Lcell,all_cell$Rcell) %>% unique(.) 
  
  ct_auc = lapply(1:length(all_cell), function(x){
    lc = all_cell[[x]]
    
    lc_enriched = selected_features %>% unlist(.) %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
      set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>%
      subset(Lcell == lc) %>% .$Interaction %>% unique(.)
    
    lc_lirics = lirics_relevant %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
      set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>%
      subset(Lcell == lc) %>% .$Interaction %>% unique(.) %>% subset(!(. %in% lc_enriched))
    
    nlc_enriched = selected_features %>% unlist(.) %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
      set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>%
      subset(Lcell != lc) %>% .$Interaction %>% unique(.)
    
    nlc_lirics = lirics_relevant %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
      set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>%
      subset(Lcell != lc) %>% .$Interaction %>% unique(.) %>% subset(!(. %in% nlc_enriched))
    
    mat = data.frame(
      "enriched" = c(length(lc_enriched), length(lc_lirics)), 
      "not_enriched" = c(length(nlc_enriched), length(nlc_lirics)),
      row.names= c("RDI", "LIRICS"),
      stringsAsFactors = F
    ) %>% as.matrix(.)
    # 4 cell matrix: 
    #         enriched  not enriched
    # RDI      S is 1    L is 0
    # LIRICS  L is 1    S is 0
    res = data.frame(p=fisher.test(mat, alternative='g')$p.value)
    row.names(res)[1] = lc
    return(res)
    
    
    
  })
  
  return(ct_auc %>% do.call(rbind,.))
  
}

fisher_Rcell <- function(selected_features){
  all_cell = selected_features %>% unlist(.) %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
    set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(GenePair=paste(L,R,sep='_')) 
  all_cell = c(all_cell$Lcell,all_cell$Rcell) %>% unique(.) 
  
  ct_auc = lapply(1:length(all_cell), function(x){
    lc = all_cell[[x]]
    
    lc_enriched = selected_features %>% unlist(.) %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
      set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>%
      subset(Rcell == lc) %>% .$Interaction %>% unique(.)
    
    lc_lirics = lirics_relevant %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
      set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>%
      subset(Rcell == lc) %>% .$Interaction %>% unique(.) %>% subset(!(. %in% lc_enriched))
    
    nlc_enriched = selected_features %>% unlist(.) %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
      set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>%
      subset(Rcell != lc) %>% .$Interaction %>% unique(.)
    
    nlc_lirics = lirics_relevant %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
      set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>%
      subset(Rcell != lc) %>% .$Interaction %>% unique(.) %>% subset(!(. %in% nlc_enriched))
    
    mat = data.frame(
      "enriched" = c(length(lc_enriched), length(lc_lirics)), 
      "not_enriched" = c(length(nlc_enriched), length(nlc_lirics)),
      row.names= c("RDI", "LIRICS"),
      stringsAsFactors = F
    ) %>% as.matrix(.)
    # 4 cell matrix: 
    #           RDI LIRICS database
    # enriched  X   X
    # !enriched X   X
    res = data.frame(p=fisher.test(mat, alternative='g')$p.value)
    row.names(res)[1] = lc
    return(res)
    
    
    
  })
  
  return(ct_auc %>% do.call(rbind,.))
  
}

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## load ICB data sets
test_input = ICB_input

## ------- Calculate Ligand Cell Type Enrichment -------
Lcell_enrichment=fisher_Lcell(feature_RDI) %>% mutate(celltype=row.names(.)) %>% arrange(p) %>% 
  mutate(rank=1:nrow(.)) %>% mutate(P=log(1/p)) %>% mutate(role='ligand')

## ------- Calculate Receptor Cell Type Enrichment -------
Rcell_enrichment=fisher_Rcell(feature_RDI) %>% mutate(celltype=row.names(.)) %>% arrange(p) %>% 
  mutate(rank=1:nrow(.)) %>% mutate(P=log(1/p)) %>% mutate(role='receptor')

## ------- Figure 3B -------
cell_enrichment = rbind(Lcell_enrichment, Rcell_enrichment)

fig3b=ggplot(data=cell_enrichment, aes(x=celltype, y=P, fill=role)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,11))+
  geom_segment(data=cell_enrichment %>% subset(role == 'ligand'), aes(x = as.numeric(factor(celltype))+0.15, y = 0, xend = as.numeric(factor(celltype))+0.15, yend = P))+
  geom_segment(data=cell_enrichment %>% subset(role == 'receptor'), aes(x = as.numeric(factor(celltype))-0.15, y = 0, xend = as.numeric(factor(celltype))-0.15, yend = P))+
  geom_point(data=cell_enrichment %>% subset(role == 'ligand'), aes(x = as.numeric(factor(celltype))+0.15, y = P), color='#ff1b6b', size=3)+
  geom_point(data=cell_enrichment %>% subset(role == 'receptor'), aes(x = as.numeric(factor(celltype))-0.15, y = P), color='#45caff', size=3)+
  theme_pubr()+
  ylab('-log(p)')+
  geom_hline(yintercept=log(1/.05), linetype="dashed", color = "black", linewidth=0.5)+ 
  theme(legend.position = 'right', axis.title.x = element_blank())+
  scale_x_continuous(breaks = as.numeric(factor(cell_enrichment$celltype)), labels = cell_enrichment$celltype)+
  theme(legend.position = c(0.85,0.7), axis.title.x = element_blank(), legend.title = element_blank())

plot(fig3b)

