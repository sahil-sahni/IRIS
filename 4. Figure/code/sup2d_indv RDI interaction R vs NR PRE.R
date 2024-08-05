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
library(ggrepel)
library(gridExtra)

# Updated from original submission 05/22/24 SS

## ------- Functions -------

fisher_lirics <- function(patient_binary, patient_group){  
  res = patient_binary[names(patient_binary) %in% subset(patient_group, response==1)$sample] #vector of all response patients
  nres = patient_binary[names(patient_binary) %in% subset(patient_group, response==0)$sample] #vector of all non-response patients
  mat = data.frame(
    "present" = c(length(subset(res, res==1)), length(subset(nres, nres==1))), 
    "absent" = c(length(subset(res, res==0)), length(subset(nres, nres==0))),
    row.names= c("response", "non response"),
    stringsAsFactors = F
  ) %>% as.matrix(.)
  # 4 cell matrix: 
  #           present absent
  # response  R is 1  R is 0
  # nonresp.  N is 1  N is 0
  res = data.frame(p=fisher.test(mat, alternative='t')$p.value, OR = fisher.test(mat)$estimate)
  
  return(res)
}
## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## load ICB data sets
model_input = ICB_input 

## load selected features
AR_LIRICS = feature_RDI
AR=AR_LIRICS %>% unlist(.)  %>% table(.) %>% as.data.frame(.) %>% set_colnames(c('interaction','freq'))

## organize input
pre=do.call(dplyr::bind_rows, model_input[[1]][c(2,3,5,6,8)]) #only pre-treatment

icb_meta = pre %>% select(sample=sample, response=response)
icb_binary = pre[,-1:-2] %>% set_rownames(pre$sample) %>% t(.) %>% .[(row.names(.) %in% AR$interaction),]

## ------- Calculate Fisher Enrichment -------
interaction_result = lapply(1:nrow(icb_binary), function(x) fisher_lirics(icb_binary[x,], icb_meta))
interaction_result = do.call(rbind,interaction_result)

AR_result = data.frame(interaction=row.names(icb_binary)) %>% cbind(., interaction_result) %>% mutate(OR=OR)

Cell_LR_df = data.frame(AR_result$interaction, str_split_fixed(AR_result$interaction, '\\_', 4)) %>%
  set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(Lcell_Rcell=paste(Lcell,Rcell,sep='_'))

AR_result = merge(AR_result, Cell_LR_df, by.x=1, by.y=1) %>% mutate(ct_FDR=NA) %>% mutate(FDR=p.adjust(p,'fdr'))

## cell-pair-specific p-value adjustment.
cell.pair.counts <- AR_result$Lcell_Rcell %>% table %>% sort(decreasing = T)       # counts of all cell pairs present
for (cc in cell.pair.counts %>% names){
  idx <- which(AR_result$Lcell_Rcell == cc)                                      # indices for the specific cell pair
  AR_result[idx, "ct_FDR"] <- AR_result[idx, "p"] %>% 
    p.adjust(method = "fdr")
}

## for source data file
AR_result = AR_result %>% mutate(group=ifelse(ct_FDR<0.2 & OR>1, 'AR','NS')) %>% arrange(ct_FDR) %>% mutate(rank=1:nrow(.))


## ------- Supp 4E -------
sup4e=ggplot(AR_result, aes(x=OR, y=log2(1/ct_FDR), color=group)) + 
  geom_point()  + 
  theme_pubr()+
  theme(legend.position = "none") +theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept=log2(1/.2), linetype="dashed", color = "grey50", linewidth=0.5)+
  geom_vline(xintercept=1, linetype="dashed", color = "grey50", linewidth=0.5)+
  scale_color_manual(values=c('#D0759F','lightgray'))+
  ylab('-log(FDR per celltype pair)')+
  ggtitle('RDI pre-treatment') + xlab('odds ratio')
plot(sup4e)
