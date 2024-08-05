library(data.table)
library(tidyverse)
#### Spatial references

## Thrane 2018 Treatment Naive:

## GEP in counts:
folder_path = '/Users/sahnis2/Documents/anti_pdl1_project/spatial data sets/Thrane_2018/ST-Melanoma-Datasets_1' 
file_names <- list.files(folder_path)

for(i in 1:length(file_names)){
  path = paste(folder_path,'/',file_names[[i]],sep="")
  
  spGEP = fread(path)
  
  spGEP = separate(spGEP, gene, into = c("GENES", "gene_identifier"), sep = " ") %>% select(-gene_identifier)
  
  # Identify unique genes
  unique_genes <- !duplicated(spGEP$GENES)
  
  # Subset the data.frame to keep only rows with unique genes
  spGEP <- spGEP[unique_genes, ] # in counts
  
  coord_df = data.frame(SpotID=colnames(spGEP)[-1])
  coord_df <- cbind(coord_df, separate(coord_df, SpotID, into = c("row", "col"), sep = "x"))
  
  sample_name <- sub("^ST_(.*?)_counts\\.tsv$", "\\1", file_names[[i]])
  
  save_path = '/Users/sahnis2/Documents/anti_pdl1_project/7. spatial analysis/cytospace/input/spatial reference_thrane_2018'

  spGEP_path = paste(save_path,'/',sample_name,'_spGEP.txt',sep='')
  spCoord_path = paste(save_path,'/',sample_name,'_spCoord.txt',sep='')
  
  write.table(spGEP, file = spGEP_path, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(coord_df, file = spCoord_path, sep = "\t", quote = FALSE, row.names = FALSE)
}
