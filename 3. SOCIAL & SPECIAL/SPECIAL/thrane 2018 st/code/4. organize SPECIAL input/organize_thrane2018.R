library(data.table)
library(tidyverse)
library(magrittr)
library(readxl)
#### Spatial references

## ------ scRNA-seq data ------
exp = fread('GSE115978_tpm.csv') %>% as.data.frame(.)
sc_gene = toupper(exp$V1)

exp = exp[,-1] %>% as.matrix(.)
rownames(exp) = sc_gene

exp = 10*((2^(exp))-1) #convert to TPM

## ------ CCI data ------

## read interaction pairs
pairs = read_xlsx('NIHMS1770717-supplement-2.xlsx', skip=1) %>% dplyr::select(ligand, receptor) # data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8983586/bin/NIHMS1770717-supplement-2.xlsx
pairs = toupper(as.matrix(pairs)) 
pairs = unique(pairs)

## filter CCI data
pairs_check = pairs %>% as.data.frame(.) %>% mutate(receptor=str_split_fixed(pairs[,'receptor'], ";", 3))  %>% as.matrix(.) 
ind1 = (pairs_check[,1] %in% sc_gene); 
ind2 = (pairs_check[,2] %in% sc_gene);
ind3 = (pairs_check[,3] %in% sc_gene)|(pairs_check[,3] == "");
ind4 = (pairs_check[,4] %in% sc_gene)|(pairs_check[,4] == "");
pairs_filtered = pairs[ind1&ind2&ind3&ind4,]

## ------ gather cell location per patient (loop) ------
folder_path <- '/Users/sahnis2/Documents/anti_pdl1_project/7. spatial analysis/cytospace/output/thrane_2018/012424' # CytoSPACE output
file_list <- list.files(folder_path)

samples = lapply(file_list, function(x){
  y = sub("([^_]+_[^_]+)_.*", "\\1", x)
  z = gsub("assigned", "", y)
})

samples = unlist(samples)

input=lapply(1:length(file_list), function(x){
  path = paste(folder_path,'/',file_list[[x]], sep='')
  
  loc = fread(path)
  
  output = list(name= samples[[x]], cci=pairs_filtered, loc=loc, exp=exp)
  return(output)
})

names(input) = samples

saveRDS(input, file='SPECIAL_input', compress = 'gzip')
