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

# Updated from original submission 05/22/24 SS

## ------- Functions -------

fisher_cellpair <- function(selected_features){
  all_celltypes = selected_features %>% unlist(.) %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
    set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(Lcell_Rcell=paste(Lcell,Rcell,sep='_'))
  all_celltypes = c(all_celltypes$Lcell, all_celltypes$Rcell) %>% unique(.) 
  
  ct_auc = lapply(1:length(all_celltypes), function(x){
    ligand_ct = all_celltypes[[x]]
    
    all_receptor_ct = all_celltypes[-x]
    
    pair_ct = lapply(1:length(all_receptor_ct), function(y){
      receptor_ct = all_receptor_ct[[y]]
      
      ct_enriched = selected_features %>% unlist(.) %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
        set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(Lcell_Rcell=paste(Lcell,Rcell,sep='_')) %>%
        subset(Lcell == ligand_ct & Rcell == receptor_ct) %>% .$Interaction %>% unique(.)
      
      ct_lirics = lirics_relevant %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
        set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(Lcell_Rcell=paste(Lcell,Rcell,sep='_')) %>%
        subset(Lcell == ligand_ct & Rcell == receptor_ct) %>% .$Interaction %>% unique(.) %>% subset(!(. %in% ct_enriched))
      
      nct_enriched = selected_features %>% unlist(.) %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
        set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(Lcell_Rcell=paste(Lcell,Rcell,sep='_')) %>%
        subset(Lcell != ligand_ct | Rcell != receptor_ct) %>% .$Interaction %>% unique(.)
      
      nct_lirics = lirics_relevant %>% data.frame(., str_split_fixed(., '\\_', 4)) %>%
        set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(Lcell_Rcell=paste(Lcell,Rcell,sep='_')) %>%
        subset(Lcell != ligand_ct | Rcell != receptor_ct) %>% .$Interaction %>% unique(.) %>% subset(!(. %in% nct_enriched))
      
      mat = data.frame(
        "enriched" = c(length(ct_enriched), length(ct_lirics)), 
        "not_enriched" = c(length(nct_enriched), length(nct_lirics)),
        row.names= c("AR", "LIRICS"),
        stringsAsFactors = F
      ) %>% as.matrix(.)
      # 4 cell matrix: 
      #         enriched  not enriched
      # RDI      S is 1    L is 0
      # LIRICS  L is 1    S is 0
      res = data.frame(p=fisher.test(mat, alternative='g')$p.value)
      
      row.names(res)[1] =paste(ligand_ct, receptor_ct, sep='_')
      return(res)
    })
    
    
    return(pair_ct %>% do.call(rbind,.))
  })
  
  return(ct_auc %>% do.call(rbind,.))
  
}
## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## load ICB data sets
test_input = ICB_input


## ------- Calculate Cell Pair Enrichment -------
cp_enrichment=fisher_cellpair(feature_RDI) %>% mutate(celltype=row.names(.)) %>% arrange(p) %>% mutate(rank=1:nrow(.)) %>% mutate(P=log2(1/p))
cp_enrichment$celltype = str_replace(cp_enrichment$celltype, '_', ' - ')
cp_enrichment$celltype[9:nrow(cp_enrichment)] <- NA
cp_enrichment$group = ifelse(cp_enrichment$p<0.05,'AR','NS')

## ------- Supp 2i -------
sup2i=ggplot(cp_enrichment, aes(x=rank, y=P,color=group)) + 
  geom_point(size = 1)+
  geom_hline(yintercept=log2(1/.05), linetype="dashed", color = "grey50", linewidth=0.5)+ 
  geom_label_repel(aes(label = celltype),
                   box.padding   = 0.35, force=25,
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   na.rm=TRUE) +
  theme_pubr()+
  theme(legend.position = "none") +
  ylab('-log2(p)')+
  scale_color_manual(values=c('#D0759F','lightgray'))

plot(sup2i)


