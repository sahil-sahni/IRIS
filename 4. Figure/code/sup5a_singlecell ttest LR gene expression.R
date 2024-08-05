library(tidyverse)
library(readxl)
library(ggrepel)
library(magrittr)
library(ggpubr)
library(data.table)
library(coin)

# New file from original submission 05/22/24 SS

## ------- Functions -------
calc_ct_LR_ttest <- function(i, meta, cci, expr){ #per cell type (i)
  ct_meta = subset(meta, celltype == i) #subset single cell for cell type
  ct_cci = subset(cci, Lcell == i | Rcell == i) # subset  cci for cell type
  ct_cci_gene = unique(c(ct_cci$L, ct_cci$R) %>% str_split(.,';')) %>% unlist(.)
  
  ct_tpm = expr %>% .[rownames(expr) %in% ct_cci_gene,c(ct_meta$Cell)] #filter single cell for cci genes and Cell barcode corresponding to cell type of interest 
  
  pout = lapply(1:nrow(ct_tpm), function(j) calc_gene_ttest(j, i, ct_tpm)) # for each gene calculate the t test 
  
  pout = do.call(rbind,pout) 
  pout = pout %>% mutate(fdr = p.adjust(pval, method='fdr'), logfdr = -log(fdr)) #perform MHC for each cell type individually
  pout = pout %>% arrange(fdr) %>% mutate(rank=1:nrow(.), significant = factor(ifelse(fdr < 0.2, 'Y','N'), levels=c('Y','N'))) 
  
  # fraction of interactions that are significant
  sig = pout %>% subset(fdr < 0.2) %>% nrow(.)
  tot = nrow(pout)
  print(paste('fraction of significantly expressed genes  in ',i, '=',sig/tot))
  
  return(pout)
}

calc_gene_ttest <- function(gene, celltype, expr){ # for each gene calculate the t test 
  
  singlegene_singlecelltype_tpm = expr[gene,] %>% t(.) %>% unlist(.)
  p =t.test(singlegene_singlecelltype_tpm, mu=0, alternative = 'greater', exact=T)$p.value # calculate p value
  if (is.nan(p)){p = 1}
  output = data.frame(celltype=celltype, gene=row.names(expr[gene,]), pval=p, mean=mean(singlegene_singlecelltype_tpm))
  
  return(output)
}

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

sc_data = fread('GSE115978_tpm.csv') %>% as.data.frame(.) # download file from: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978%5Ftpm%2Ecsv%2Egz

# load PCG 
pcg = pcg #protein coding genes

# load livnat meta files 
sc_celltype = livnat_celltype
sc_meta = livnat_meta

## 1. convert livnat log10(TPM) to TPM
rownames(sc_data) = sc_data$V1
sc_data = sc_data[rownames(sc_data) %in% pcg$V1,-1] # subset for protein coding genes
sc_data = 10*((2^(sc_data))-1) # convert to TPM with PCG only

## 2. reorganize livnat meta
sc_celltype$celltype = plyr::mapvalues(sc_celltype$celltype, from=c('B.cell','Endo.','T.CD8','T.CD4'), to=c('Bcell','Endo','TCD8','TCD4'))
sc_meta = sc_meta %>% subset(Patient != 'Mel04-3') %>% # remove the one post-ICB responder
  mutate(Group = ifelse(Treatment == 'None', 'Naïve','Resistant')) %>% select(Cell, Group) 
sc_meta_merge = merge(sc_meta, sc_celltype, by.x='Cell', by.y='cell')

## 3.  Load Interactions
RDI = feature_RDI
RDI = RDI[c(2,3,5,6,8)] %>% flatten(.) %>% plyr::ldply(., cbind) %>% set_colnames(c('cohort','interaction')) %>% data.frame(.,str_split_fixed(.$interaction, '\\_', 4)) %>%
  set_colnames(c('cohorts','interaction','Lcell','Rcell','L','R')) %>% mutate(GenePair=paste(L,R,sep='_')) 
RDI = RDI %>% select(2:7) %>% distinct(.)

## ------- Calculate one sided t test for each gene for each cell type -------
ttest_res = lapply(unique(sc_meta_merge$celltype), function(i) calc_ct_LR_ttest(i, sc_meta_merge, RDI, sc_data)) #per cell type
ttest_res = do.call(rbind, ttest_res)

# calculate the fractions of genes that are signficant
ttest_frac = ttest_res %>% mutate(significant = ifelse(significant == "Y", 1, 0)) 
frac_significant <- ttest_frac %>%
  group_by(celltype) %>%
  summarize(fraction_significant = mean(significant)*100)
  
## ------- Supp 5a -------
sup5a = ggplot(ttest_res, aes(rank,logfdr))+
  geom_point(aes(color=significant),size=2)+
  theme_pubr()+
  scale_color_manual(values=c('#D0759F','lightgray'))+
  ylab('-log(FDR)')+
  geom_hline(yintercept=log(1/.2), linetype="dashed", color = "grey50", linewidth=0.5)+
  theme(legend.position = c('none'))+
  facet_wrap(~celltype,nrow=2,ncol=5,scales = 'free')+
  theme(strip.text.x = element_text(size = 15))


sup5a= sup5a + geom_text(data = frac_significant, aes(x = Inf, y = Inf, label = paste0("% Significant: ", round(fraction_significant,1))), 
                   hjust = 1, vjust = 2, size = 5)

plot(sup5a)

