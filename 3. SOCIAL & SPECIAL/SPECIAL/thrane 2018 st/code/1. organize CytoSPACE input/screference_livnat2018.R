library(data.table)
library(tidyverse)
#### Single-cell references

## Livnat 2018 Treatment Naive:

## cell information
cl = readRDS('Livnat_updated_celltypes_table.rds')
cl = cl %>% subset(celltype %in% c('Mal','T.CD8','B.cell','T.CD4','Macrophage','Endo.','CAF','NK','skinDC','pDC'))

## meta information:
meta = fread('livnat_meta_tisch.txt') # download from 
filtered_cell = meta %>% subset(Treatment == 'None' & Cell %in% cl$cell) %>% .$Cell
filtered_cell = meta$Cell

## scRNA-seq gene expression file (raw counts):
scgep = fread('GSE115978_counts.csv')
colnames(scgep)[[1]] = 'GENES'

# subset for only treatment naive cells
scgep = scgep %>% as.data.frame(.) %>% select(c('GENES', filtered_cell))
cl = cl %>% subset(cell %in% colnames(scgep)) %>% select(`Cell IDs`=1, CellType=2) %>% mutate(CellType = plyr::mapvalues(CellType, from=c('B.cell','Endo.','T.CD8','T.CD4'), to=c('Bcell','Endo','TCD8','TCD4')))

#save file
write.table(scgep, file = "livnat18_scGEP_naive.txt", sep = "\t", quote = FALSE, row.names=F)
write.table(cl, file = "livnat18_scCellLabels_naive.txt", sep = "\t", quote = FALSE, row.names = F)


