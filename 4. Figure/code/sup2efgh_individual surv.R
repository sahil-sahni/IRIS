library(tidyverse)
library(magrittr)
library(plyr)
library(readxl)
library(ModelMetrics)
library(aod)
library(cutpointr)
library(survival)
library(ggpubr)
library(survminer)
###

# New file from original submission 05/23/24 SS

## ------- Functions -------
calc_HR_cci = function(z, cci_list, cci_df){
  
  cci = cci_df[,cci_list[[z]]] # subset for the binary active and inactive activity from df
  surv_input = data.frame(sample=cci_df$sample, status=as.numeric(cci_df[,3]), days=cci_df[,4], cci=cci)
  surv_input = na.omit(surv_input)
  cox.out = coxph(Surv(days,status) ~ cci, data=surv_input)
  aa  = summary(cox.out)        
  output = aa$coefficients["cci",c(1,5)]
  
  return(output)
}

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

feature_RDI_indv = unlist(feature_RDI) %>% unique(.)
feature_RUI_indv = unlist(feature_RUI) %>% unique(.)

## ------- ################ -------
## ------- ICB OVERALL SURVIVAL (MERGE) RDI -------
## ------- ################ -------

## ------- Input -------

## load external ICB data sets
merged_input_OS = ICB_OS_input[[1]]
merged_input_OS = do.call(rbind.fill, merged_input_OS[c(2,3,5,6)])

## ------- Calculate Hazard Ratio -------
OS_HR_RDI_res = lapply(1:length(feature_RDI_indv), function(z) calc_HR_cci(z, feature_RDI_indv, merged_input_OS))
OS_HR_RDI_res = do.call(rbind, OS_HR_RDI_res) 
row.names(OS_HR_RDI_res) = feature_RDI_indv

OS_HR_RDI_res = OS_HR_RDI_res %>% as.data.frame(.) %>% set_colnames(c('HR','p')) %>% mutate(interaction = rownames(.))

# format data frame for multiple hypothesis testing
Cell_LR_df = data.frame(OS_HR_RDI_res$interaction, str_split_fixed(OS_HR_RDI_res$interaction, '\\_', 4)) %>%
  set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(Lcell_Rcell=paste(Lcell,Rcell,sep='_'))

OS_HR_RDI_res = merge(OS_HR_RDI_res, Cell_LR_df, by.x=3, by.y=1) %>% mutate(ct_FDR=NA) %>% mutate(FDR=p.adjust(p,'fdr'))

# cell-pair-specific p-value adjustment.
cell.pair.counts <- OS_HR_RDI_res$Lcell_Rcell %>% table %>% sort(decreasing = T)       # counts of all cell pairs present
for (cc in cell.pair.counts %>% names){
  idx <- which(OS_HR_RDI_res$Lcell_Rcell == cc)                                      # indices for the specific cell pair
  OS_HR_RDI_res[idx, "ct_FDR"] <- OS_HR_RDI_res[idx, "p"] %>% 
    p.adjust(method = "fdr")
}

OS_HR_RDI_res = OS_HR_RDI_res %>% mutate(group=ifelse(ct_FDR<0.2 & HR < 0, 'AR','NS'))

## ------- for Supp Table 1 -------

st1 = OS_HR_RDI_res %>% dplyr::select(2,9)
row.names(st1) = OS_HR_RDI_res$interaction

## ------- Supp Figure 2G -------

sup2g=ggplot(OS_HR_RDI_res, aes(x=HR, y=log2(1/ct_FDR), color=group)) + 
  geom_point()  + 
  theme_pubr()+
  theme(legend.position = "none") +theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept=log2(1/.2), linetype="dashed", color = "grey50", linewidth=0.5)+
  geom_vline(xintercept=0, linetype="dashed", color = "grey50", linewidth=0.5)+
  scale_color_manual(values=c('#D0759F','lightgray'))+
  ylab('-log(FDR per cell type pair)')+
  ggtitle('RDI overall survival') + xlab('hazard ratio')
plot(sup2g)

## ------- ################ -------
## ------- ICB PROGRESSION FREE SURVIVAL (MERGE) RDI -------
## ------- ################ -------

## ------- Input -------

## load external ICB data sets
merged_input_PFS = ICB_PFS_input[[1]]
merged_input_PFS = do.call(rbind.fill, merged_input_PFS[c(2,3,5,6)])

## ------- Calculate Hazard Ratio -------
PFS_HR_RDI_res = lapply(1:length(feature_RDI_indv), function(z) calc_HR_cci(z, feature_RDI_indv, merged_input_PFS))
PFS_HR_RDI_res = do.call(rbind, PFS_HR_RDI_res) 
row.names(PFS_HR_RDI_res) = feature_RDI_indv

PFS_HR_RDI_res = PFS_HR_RDI_res %>% as.data.frame(.) %>% set_colnames(c('HR','p')) %>% mutate(interaction = rownames(.))

# format data frame for multiple hypothesis testing
Cell_LR_df = data.frame(PFS_HR_RDI_res$interaction, str_split_fixed(PFS_HR_RDI_res$interaction, '\\_', 4)) %>%
  set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(Lcell_Rcell=paste(Lcell,Rcell,sep='_'))

PFS_HR_RDI_res = merge(PFS_HR_RDI_res, Cell_LR_df, by.x=3, by.y=1) %>% mutate(ct_FDR=NA) %>% mutate(FDR=p.adjust(p,'fdr'))

## cell-pair-specific p-value adjustment.
cell.pair.counts <- PFS_HR_RDI_res$Lcell_Rcell %>% table %>% sort(decreasing = T)       # counts of all cell pairs present
for (cc in cell.pair.counts %>% names){
  idx <- which(PFS_HR_RDI_res$Lcell_Rcell == cc)                                      # indices for the specific cell pair
  PFS_HR_RDI_res[idx, "ct_FDR"] <- PFS_HR_RDI_res[idx, "p"] %>% 
    p.adjust(method = "fdr")
}

PFS_HR_RDI_res = PFS_HR_RDI_res %>% mutate(group=ifelse(ct_FDR<0.2 & HR < 0, 'AR','NS'))

## ------- Supp Table 1 -------

st1 = PFS_HR_RDI_res %>% dplyr::select(2,9)
row.names(st1) = PFS_HR_RDI_res$interaction

## ------- Supp Figure 2H -------

sup2h=ggplot(PFS_HR_RDI_res, aes(x=HR, y=log2(1/ct_FDR), color=group)) + 
  geom_point()  + 
  theme_pubr()+
  theme(legend.position = "none") +theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept=log2(1/.2), linetype="dashed", color = "grey50", linewidth=0.5)+
  geom_vline(xintercept=0, linetype="dashed", color = "grey50", linewidth=0.5)+
  scale_color_manual(values=c('#D0759F','lightgray'))+
  ylab('-log(FDR per cell type pair)')+
  ggtitle('RDI progression free survival') + xlab('hazard ratio')
plot(sup2h)


###################
###################
###################
###################

## ------- ################ -------
## ------- ICB OVERALL SURVIVAL (MERGE) RUI -------
## ------- ################ -------

## ------- Input -------

## load external ICB data sets
merged_input_OS = ICB_OS_input[[1]]
merged_input_OS = do.call(rbind.fill, merged_input_OS[c(2,3,5,6)])

## ------- Calculate Hazard Ratio -------
OS_HR_RUI_res = lapply(1:length(feature_RUI_indv), function(z) calc_HR_cci(z, feature_RUI_indv, merged_input_OS))
OS_HR_RUI_res = do.call(rbind, OS_HR_RUI_res) 
row.names(OS_HR_RUI_res) = feature_RUI_indv

OS_HR_RUI_res = OS_HR_RUI_res %>% as.data.frame(.) %>% set_colnames(c('HR','p')) %>% mutate(interaction = rownames(.))

# format data frame for multiple hypothesis testing
Cell_LR_df = data.frame(OS_HR_RUI_res$interaction, str_split_fixed(OS_HR_RUI_res$interaction, '\\_', 4)) %>%
  set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(Lcell_Rcell=paste(Lcell,Rcell,sep='_'))

OS_HR_RUI_res = merge(OS_HR_RUI_res, Cell_LR_df, by.x=3, by.y=1) %>% mutate(ct_FDR=NA) %>% mutate(FDR=p.adjust(p,'fdr'))

## cell-pair-specific p-value adjustment.
cell.pair.counts <- OS_HR_RUI_res$Lcell_Rcell %>% table %>% sort(decreasing = T)       # counts of all cell pairs present
for (cc in cell.pair.counts %>% names){
  idx <- which(OS_HR_RUI_res$Lcell_Rcell == cc)                                      # indices for the specific cell pair
  OS_HR_RUI_res[idx, "ct_FDR"] <- OS_HR_RUI_res[idx, "p"] %>% 
    p.adjust(method = "fdr")
}

OS_HR_RUI_res = OS_HR_RUI_res %>% mutate(group=ifelse(ct_FDR<0.2 & HR > 0, 'NR','NS'))


## ------- Supp Figure 2E -------

sup2e=ggplot(OS_HR_RUI_res, aes(x=HR, y=log2(1/ct_FDR), color=group)) + 
  geom_point()  + 
  theme_pubr()+
  theme(legend.position = "none") +theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept=log2(1/.2), linetype="dashed", color = "grey50", linewidth=0.5)+
  geom_vline(xintercept=0, linetype="dashed", color = "grey50", linewidth=0.5)+
  scale_color_manual(values=c('#88BECD','lightgray'))+
  ylab('-log(FDR per cell type pair)')+
  ggtitle('RUI overall survival') + xlab('hazard ratio')
plot(sup2e)



## ------- ################ -------
## ------- ICB PROGRESSION FREE SURVIVAL (MERGE) RUI -------
## ------- ################ -------

## ------- Input -------

## load external ICB data sets
merged_input_PFS = ICB_PFS_input[[1]]
merged_input_PFS = do.call(rbind.fill, merged_input_PFS[c(2,3,5,6)])

## ------- Calculate Hazard Ratio -------
PFS_HR_RUI_res = lapply(1:length(feature_RUI_indv), function(z) calc_HR_cci(z, feature_RUI_indv, merged_input_PFS))
PFS_HR_RUI_res = do.call(rbind, PFS_HR_RUI_res) 
row.names(PFS_HR_RUI_res) = feature_RUI_indv

PFS_HR_RUI_res = PFS_HR_RUI_res %>% as.data.frame(.) %>% set_colnames(c('HR','p')) %>% mutate(interaction = rownames(.))

# format data frame for multiple hypothesis testing
Cell_LR_df = data.frame(PFS_HR_RUI_res$interaction, str_split_fixed(PFS_HR_RUI_res$interaction, '\\_', 4)) %>%
  set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(Lcell_Rcell=paste(Lcell,Rcell,sep='_'))

PFS_HR_RUI_res = merge(PFS_HR_RUI_res, Cell_LR_df, by.x=3, by.y=1) %>% mutate(ct_FDR=NA) %>% mutate(FDR=p.adjust(p,'fdr'))

## cell-pair-specific p-value adjustment.
cell.pair.counts <- PFS_HR_RUI_res$Lcell_Rcell %>% table %>% sort(decreasing = T)       # counts of all cell pairs present
for (cc in cell.pair.counts %>% names){
  idx <- which(PFS_HR_RUI_res$Lcell_Rcell == cc)                                      # indices for the specific cell pair
  PFS_HR_RUI_res[idx, "ct_FDR"] <- PFS_HR_RUI_res[idx, "p"] %>% 
    p.adjust(method = "fdr")
}

PFS_HR_RUI_res = PFS_HR_RUI_res %>% mutate(group=ifelse(ct_FDR<0.2 & HR > 0, 'NR','NS'))

## ------- Supp Figure 2F -------

sup2f=ggplot(PFS_HR_RUI_res, aes(x=HR, y=log2(1/ct_FDR), color=group)) + 
  geom_point()  + 
  theme_pubr()+
  theme(legend.position = "none") +theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept=log2(1/.2), linetype="dashed", color = "grey50", linewidth=0.5)+
  geom_vline(xintercept=0, linetype="dashed", color = "grey50", linewidth=0.5)+
  scale_color_manual(values=c('#88BECD','lightgray'))+
  ylab('-log(FDR per cell type pair)')+
  ggtitle('RUI progression free survival') + xlab('hazard ratio')
plot(sup2f)



