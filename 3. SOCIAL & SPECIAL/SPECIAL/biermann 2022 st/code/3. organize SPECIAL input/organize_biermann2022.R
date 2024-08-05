library(data.table)
library(tidyverse)
library(magrittr)
library(readxl)
library(Matrix)
#### Spatial references

## ------ gather cell location per patient (loop) ------
folder_path <- '/Users/sahnis2/Documents/anti_pdl1_project/7. spatial analysis/cytospace/output/izar_2022/030224_v2' #CytoSPACE output
file_list <- list.files(folder_path)

scfolder_path = '/Users/sahnis2/Documents/anti_pdl1_project/spatial data sets/Biermann_2022/GSE185386_snseq_match' #expression data of matched snRNA-seq and spatial
scfile_list <- list.files(scfolder_path)


samples = lapply(file_list, function(x){
  y <- sub("^([A-Za-z]+\\d+\\.\\d)_.*", "\\1", x)
  z <- gsub("assigned", "", y)
})

samples = unlist(samples) %>% unique(.)

input=lapply(1:length(samples), function(x){
  ## ------ spatial data ------
  ## first cell type path
  path = paste(folder_path,'/',samples[[x]],'_firstctassigned_locations.csv', sep='')
  firstloc = fread(path)
  
  ## second cell type path
  path = paste(folder_path,'/',samples[[x]],'_secondctassigned_locations.csv', sep='')
  secondloc = fread(path)
  
  ## combine cell type location
  loc = rbind(firstloc, secondloc)
  
  ## ------ scRNA-seq data ------
  ## extract matched single cell data
  i = gsub("\\.\\d+", "", samples[[x]])
  scname = plyr::mapvalues(i, to=c("MPM01", "MPM02", "MPM04", "MPM05", "MPM06", "MPM07", "MPM08", "MPM09", "MPM10", "MPM11"),
                      from=c("ECM01", "ECM02", "ECM04", "ECM05", "ECM06", "ECM07", "ECM08", "ECM09", "ECM10", "ECM11"))
  
  matchedsc_filepath <- paste(scfolder_path,'/',scfile_list[grep(scname, scfile_list)],sep="")
  exp=fread(matchedsc_filepath)
  
  ## reformat sc data
  sc_gene = toupper(exp$V1)
  
  exp = exp[,-1] %>% as.matrix(.)
  rownames(exp) = sc_gene
  
  ## convert to TPM
  exp = t(10^6 * (t(exp) / Matrix::colSums(exp)))
  print(colSums(exp))
  
  ## ------ CCI data ------
  ## read interaction pairs
  pairs = read_xlsx('NIHMS1770717-supplement-2.xlsx', skip=1) %>% dplyr::select(ligand, receptor) #data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8983586/bin/NIHMS1770717-supplement-2.xlsx
  pairs = toupper(as.matrix(pairs))
  pairs = unique(pairs)
  
  ## filter CCI data
  pairs_check = pairs %>% as.data.frame(.) %>% mutate(receptor=str_split_fixed(pairs[,'receptor'], ";", 3))  %>% as.matrix(.) 
  ind1 = (pairs_check[,1] %in% sc_gene); 
  ind2 = (pairs_check[,2] %in% sc_gene);
  ind3 = (pairs_check[,3] %in% sc_gene)|(pairs_check[,3] == "");
  ind4 = (pairs_check[,4] %in% sc_gene)|(pairs_check[,4] == "");
  pairs_filtered = pairs[ind1&ind2&ind3&ind4,]

  output = list(name= samples[[x]], cci=pairs_filtered, loc=loc, exp=exp)
  return(output)
})

names(input) = samples

saveRDS(input, file='SPECIAL_input.RDS')
