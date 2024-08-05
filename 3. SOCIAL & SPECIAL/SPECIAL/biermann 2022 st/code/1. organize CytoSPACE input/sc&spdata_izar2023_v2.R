library(data.table)
library(tidyverse)
library(singleCellNet)
library(magrittr)
library(spacexr)
library(Matrix)

save_path = '/Users/sahnis2/Documents/anti_pdl1_project/7. spatial analysis/cytospace/input/spatial reference_izar_2022_v2'

#### Generate inputs for Izar et al. 2023 Slide-seq V2 ####

### 2. prepare "queried data" Izar et al. 2023 data ###

# sc meta data

scmeta = fread('/Users/sahnis2/Documents/anti_pdl1_project/spatial data sets/Biermann_2022/GSE200218_sc_sn_metadata.csv') %>% as.data.frame() %>% select(V1,ID,cell_type_int)
orgct = c('Plasma cells', 'Tumor cells', 'Doublets', 'Low-quality cells', 'Cycling cells' , 'Dendritic cells', 'Microglia', 'MDM', 'Monocytes', 'Tregs', 'CD8+ T cells', 'B cells', 'NK cells', 'CD4+ T cells', 'CNS cells', 'Stromal cells', 'Endothelial cells', 'Contamination', 'Mast cells', 'Epithelial cells')
react = c('Bcell', 'Mal', 'OMIT', 'OMIT', 'OMIT', 'skinDC', 'Macrophage', 'Macrophage', 'Additional-Immune', 'Additional-Immune', 'TCD8', 'Bcell', 'NK', 'TCD4', 'CNS', 'CAF', 'Endo', 'OMIT', 'Additional-Immune', 'Epithelial')
scmeta$cell_type_reannotate = plyr::mapvalues(scmeta$cell_type_int, from=orgct, to=react)

# sc 10X (all) file path:

scfolder_path = '/Users/sahnis2/Documents/anti_pdl1_project/spatial data sets/Biermann_2022/GSE185386_snseq_all'
scfile_names <- list.files(scfolder_path)

# sc 10X (matched) file path:

scmfolder_path = '/Users/sahnis2/Documents/anti_pdl1_project/spatial data sets/Biermann_2022/GSE185386_snseq_match'
scmfile_names <- list.files(scmfolder_path)[c(-11,-1:-3)] #remove 1
scm_samples = str_extract(scmfile_names, "(MBM|MPM)\\d+")

# slide-seqV2 file path:

spfolder_path = '/Users/sahnis2/Documents/anti_pdl1_project/spatial data sets/Biermann_2022/GSE185386_slideseq'
spfile_names <- list.files(spfolder_path)
sp_samples <- str_remove_all(spfile_names, "_slide_raw_counts.csv.gz") %>% str_remove_all(., "_slide_spatial_info.csv.gz") %>% unique(.)

### 3. Merge snRNA-seq data for reannotation purposes and robust cell fraction estimation

sc10X_merge = NULL
for(i in scfile_names){
  print(i)
  scpath = paste(scfolder_path,'/',i,sep="")
  sc10X_i = fread(scpath) %>% as.data.frame(.) #load file
  
  genename = sc10X_i$V1 #keep gene name
  sc10X_i = sc10X_i[,-1] %>% as.matrix(.) %>% Matrix::Matrix(.,sparse=T) #convert to 10X counts


  rownames(sc10X_i)=genename #replace rowname with kept genes
  
  sc10X_merge=cbind(sc10X_merge,sc10X_i)
}

### 4. Filter cell type

cl =  scmeta  %>%
  subset(cell_type_reannotate != 'OMIT') %>% # remove omitted cell type
  subset(cell_type_reannotate != 'Additional-Immune') %>% # remove Monocytes
  select(`Cell IDs`=1, CellType=4, sample=2) 

filter_ct = names(table(cl$CellType)[which(table(cl$CellType) > 25)]) #filter for cell types with > 25
cl = cl %>% subset(CellType %in% filter_ct) %>% subset(`Cell IDs` %in% colnames(sc10X_merge))

# RCTD cell input 
counts_sc = sc10X_merge[,cl$`Cell IDs`]
cell_types = factor(cl$CellType)
names(cell_types) = cl$`Cell IDs`


### 3. loop to create input for cytospace and create input for spatial analysis ###

for (i in scm_samples){
  
  # get single cell data
  msc_label = cl %>% subset(sample == i) %>% select(-sample)
  msc_10X = sc10X_merge[, msc_label$`Cell IDs`] %>% as.matrix(.) %>% as.data.frame(.)
  msc_10X_gene = row.names(msc_10X)
  msc_10X = cbind(msc_10X_gene, msc_10X)
  colnames(msc_10X)[[1]] = 'GENES'
  
  # do a mapvalues for spatial spots
  i = plyr::mapvalues(i, from=c("MPM01", "MPM02", "MPM04", "MPM05", "MPM06", "MPM07", "MPM08", "MPM09", "MPM10", "MPM11"),
                to=c("ECM01", "ECM02", "ECM04", "ECM05", "ECM06", "ECM07", "ECM08", "ECM09", "ECM10", "ECM11"))
                          
  # get spatial data of matched single cell sample
  sp_sample_intersect <- sp_samples[grep(i, sp_samples)]

  for (j in 1:length(sp_sample_intersect)){
    
    # rename sample
    sample_name = paste(i,j,sep='.')
    
    
    # load queried spatial counts matrix
    path = paste(spfolder_path,'/',sp_sample_intersect[[j]],'_slide_raw_counts.csv.gz',sep='')
    spGEP = fread(path)
    colnames(spGEP) <- gsub("-", "_", colnames(spGEP))
    colnames(spGEP) = paste(sample_name, colnames(spGEP), sep='_')
    colnames(spGEP)[[1]] = 'GENES'
    
    
    # load queried spatial coordinate information
    path = paste(spfolder_path,'/',sp_sample_intersect[[j]],'_slide_spatial_info.csv.gz',sep='')
    spdf = fread(path)
    spdf = data.frame(SpotID=spdf$barcode, row=spdf$xcoord, col=spdf$ycoord) %>% mutate(SpotID = gsub("-", "_", SpotID)) %>% mutate(SpotID = paste(sample_name, SpotID, sep='_'))
    
    # spatial input
    coords = data.frame(x=spdf$row, y=spdf$col, row.names = spdf$SpotID)
    counts_sp = spGEP %>% as.data.frame(.) %>% set_rownames(.$GENES) %>% .[,-1] %>% as.matrix(.)
  
    #1. Create the Reference object
    reference <- spacexr::Reference(counts_sc, cell_types)
    
    #2.  Create SpatialRNA object
    puck <- SpatialRNA(coords, counts_sp)
    
    #3. run RCTD
    set.seed(20892)
    myRCTD <- create.RCTD(puck, reference, max_cores = 8)
    myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
    
    #4. get results
    results <- myRCTD@results
    results =results$results_df
    
    #5. get results
    singlet = results %>% subset(spot_class %in% c('singlet','doublet_uncertain', 'doublet_certain')) %>% mutate(SpotID = rownames(.)) %>% select(SpotID=SpotID, CellType = first_type)
    doublet = results %>% subset(spot_class == 'doublet_certain') %>% mutate(SpotID = rownames(.)) %>% select(SpotID=SpotID, CellType = second_type)
    
    #6. remove spatial spots corresponding to non present cell types!
    singlet = singlet %>% subset(CellType %in% msc_label$CellType)
    doublet = doublet %>% subset(CellType %in% msc_label$CellType)
    
    spCTloc = rbind(singlet, doublet)
    

    
    # SAVE DATA
    spGEP_path = paste(save_path,'/',sample_name,'_spGEP.txt',sep='')
    spCoord_path = paste(save_path,'/',sample_name,'_spCoord.txt',sep='')
    spCTloc_path = paste(save_path,'/',sample_name,'_spCTloc.txt',sep='')
    label_path = paste(save_path,'/',sample_name,'_scLabel.txt',sep='')
    scGEP_path = paste(save_path,'/',sample_name,'_scGEP.txt',sep='')

    # Togeather 
    fwrite(as.data.frame(spGEP) %>% .[,c('GENES', spCTloc$SpotID)], file = spGEP_path, sep = "\t", quote = FALSE, row.names = FALSE)
    fwrite(spdf %>% subset(SpotID %in% spCTloc$SpotID), file = spCoord_path, sep = "\t", quote = FALSE, row.names = FALSE)
    fwrite(spCTloc, file = spCTloc_path, sep = "\t", quote = FALSE, row.names = FALSE)
    fwrite(msc_label, file = label_path, sep = "\t", quote = FALSE, row.names = FALSE)
    fwrite(msc_10X, file = scGEP_path, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Singlet Mode
    spGEP_path = paste(save_path,'/',sample_name,'_spGEP_firstct.txt',sep='')
    spCoord_path = paste(save_path,'/',sample_name,'_spCoord_firstct.txt',sep='')
    spCTloc_path = paste(save_path,'/',sample_name,'_spCTloc_firstct.txt',sep='')
    
    fwrite(as.data.frame(spGEP) %>% .[,c('GENES', singlet$SpotID)], file = spGEP_path, sep = "\t", quote = FALSE, row.names = FALSE)
    fwrite(spdf %>% subset(SpotID %in% singlet$SpotID), file = spCoord_path, sep = "\t", quote = FALSE, row.names = FALSE)
    fwrite(singlet, file = spCTloc_path, sep = "\t", quote = FALSE, row.names = FALSE)

    # Doublet Mode
    spGEP_path = paste(save_path,'/',sample_name,'_spGEP_secondct.txt',sep='')
    spCoord_path = paste(save_path,'/',sample_name,'_spCoord_secondct.txt',sep='')
    spCTloc_path = paste(save_path,'/',sample_name,'_spCTloc_secondct.txt',sep='')
    
    fwrite(as.data.frame(spGEP) %>% .[,c('GENES', doublet$SpotID)], file = spGEP_path, sep = "\t", quote = FALSE, row.names = FALSE)
    fwrite(spdf %>% subset(SpotID %in% doublet$SpotID), file = spCoord_path, sep = "\t", quote = FALSE, row.names = FALSE)
    fwrite(doublet, file = spCTloc_path, sep = "\t", quote = FALSE, row.names = FALSE)
    
  }
}
