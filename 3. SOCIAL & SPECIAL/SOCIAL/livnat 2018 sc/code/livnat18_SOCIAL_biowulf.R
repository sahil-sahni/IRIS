library(readxl)
library(tidyverse)
library(stringr)
library(reshape2)
library(magrittr)
library(plyr)
library(rslurm)
library(dplyr)
library(parallel)
library(data.table)

# New file from original submission 05/22/24 SS

## ------- Functions -------
## downsamplepp(): function used to randomy down sample 40% (default) of cells from each cell type for each group (i.e. post-ICI or naive)
##                 and generates two (depends on # of groups) pseudopatients (i.e. 1 post-ICI and 1 naive)

  ## input:
    ## i: nth pseudopatient (vector)
    ## ct_info: scRNA-seq cell type information (data.frame) with columns:
        ## Cell: single-cell barcode
        ## group: i.e. corresponding to what treatment group the cell corresponds to (i.e. treatment naive vs. resistant)
        ## cell_type: the cell assigned cell type
    ## prop: proportion of cells to downsample of each cell type (default is 0.4)

  ## output:
    ## output: modified scRNA-seq data.frame with columns:
      ## Cell
      ## cell_type
      ## samples: name of assigned pseudosample

downsamplepp <- function(i, ct_info, prop=0.4){
  output = ct_info %>% group_by(group, cell_type) %>% slice_sample(prop=prop) %>% ungroup(.) %>% mutate(samples=paste(group,i,sep="")) %>% dplyr::select(Cell, cell_type, samples)
  return(output)
}

## ------- Input -------
setwd("/data/sahnis2/SOCIAL/test_livnat") #CHANGE

# path to SOCIAL function (for RSLURM)
path_social = 'SOCIAL_SPECIAL.R' #CHANGE
source(path_social) #CHANGE

# data from livnat jerby arnon Cell 2018
sc_exp = read.csv("GSE115978_tpm.csv", row.names=1) %>% as.matrix(.) # download from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978%5Ftpm%2Ecsv%2Egz 
patient_meta = read_excel("NIHMS1508390-supplement-8.xlsx",sheet = "TableS1A_scCohort") # download from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6410377/bin/NIHMS1508390-supplement-8.xlsx
ct_info = readRDS("Livnat_updated_celltypes_table.rds") 
  
# data from kun wang Cancer Discovery 2022
lirics.db = read_excel("NIHMS1770717-supplement-2.xlsx",sheet = "ligand_receptor_interactions", skip=1) %>% dplyr::select(ligand, receptor) # download from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8983586/bin/NIHMS1770717-supplement-2.xlsx

## ------- Initialize environment -------
## set seed
set.seed(20892) #NIH zip code

## parallelize code
ncpus <- parallelly::availableCores()
print(ncpus)
options(mc.cores = ncpus) # set a global option for parallel packages

## ------- reformat input for SOCIAL (Change as needed) -------

#### scRNA-seq data ####
# reformat single-cell expression profile
row.names(sc_exp) = toupper(row.names(sc_exp))
sc_gene = row.names(sc_exp)

# correct scRNA-seq to TPM; Jerby-Arnon et al. published values in log2(TPM/10 + 1)
sc_exp = 10*((2^(sc_exp))-1)

#### patient-specific meta information ####
## format: data.frame with col: patient; treatment (Untreated (NAV) or Post-ICI resistant (RES))

patient_meta = patient_meta %>% set_colnames(., .[2,]) %>% .[-1:-2,-1] %>% 
  dplyr::select(patient="Sample", treatment="Treatment group") %>% 
  mutate(patient=toupper(patient), 
         treatment= plyr::mapvalues(treatment, from=c("Untreated", "Post-ICI (resistant)"), to=c('NAV','RES'))) #simplify label as NAV and RES

#### single-cell cell type specific information ####
## format: data.frame with col: Cell (barcode for sc); cell_type (updated cell type information); patient

ct_info = ct_info %>% select(Cell=1, cell_type=2, patients=3) %>% 
  magrittr::set_rownames(NULL) %>% 
  mutate(Cell=as.character(Cell), 
         cell_type=as.character(cell_type), 
         patients=toupper(as.character(patients))) # make basic ct_map outline

## make single-cell cell type compatible with deconvolved cell types
sc_ct= c('Mal', 'skinDC', 'Endo.', 'Macrophage', 'T.CD4', 'CAF', 'T.CD8', 'NK', 'B.cell', 'pDC')
deconv_ct= c('Mal', 'skinDC', 'Endo', 'Macrophage', 'TCD4', 'CAF', 'TCD8', 'NK', 'Bcell', 'pDC')

ct_info = ct_info %>% mutate(cell_type=mapvalues(cell_type,from=sc_ct, to=deconv_ct))

## assign Cell ("barcode") information with treatment group ('NAV' is naive, 'RES' is post-ICI resistant)
ct_info$group = mapvalues(ct_info$patients, from=patient_meta$patient, to=patient_meta$treatment) 
ct_info = ct_info %>% subset(group != 'OR') # remove one responder patient

## ------- create pseudopatients from single-cell meta information -------

# number of pseudopatients
npp = 200
prop = 0.4 

# generate pseudopatients
ct_map = lapply(1:npp, function(i) downsamplepp(i, ct_info))
ct_map = do.call(rbind, ct_map) %>% as.data.frame(.) #creates a ct_map with all pseudo patients

## ------- SOCIAL -------

#### STEP 1: ligand-receptor genes ####
# subset genes within ligand-receptor complexes that are shared with LIRICS database

pairs = SOCIAL.query_LRdb(lirics.db, sc_exp)

#### STEP 2&3: calculate interaction score per patient (using rslurm) ####

interaction_score = SOCIAL.cis_rslurm(sc_exp, ct_map, pairs, n_iterations = 100, path_social)
saveRDS(interaction_score, "livnat2018_TPM_400pseudopat_40ds_SOCIAL.rds")
