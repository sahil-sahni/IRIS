library(readxl)
library(data.table)

## SET PATH TO INPUT FOLDER HERE 
PATH = '~/Documents/NCI/anti_pdl1_project/documents/zenodo/nature communications revision 2/IRIS data' #CHANGE
PATH_code = '~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2' #CHANGE

## code path ##

# setwd for figures
setwd(paste(PATH_code,'/4. Figure/Output',sep=''))

## load major code
source(paste(PATH_code,'/4. Figure/Code/calculate_score.R',sep='')) #scoring function
source(paste(PATH_code,'/3. SOCIAL & SPECIAL/SOCIAL & SPECIAL Function/SOCIAL_SPECIAL.R',sep='')) #social special function

## LIRICS ## 

## load bulk ICB data sets
ICB_input = loadRData(paste(PATH,'/LIRICS/organized data/ICB/IRIS_LIRICS_final_ICB_cohort_input_response.Rdata',sep='')) #new meta
ICB_OS_input = loadRData(paste(PATH,'/LIRICS/organized data/ICB/IRIS_LIRICS_final_ICB_cohort_input_OS.Rdata',sep='')) #new meta
ICB_PFS_input = loadRData(paste(PATH,'/LIRICS/organized data/ICB/IRIS_LIRICS_final_ICB_cohort_input_PFS.Rdata',sep='')) #new meta

## load TCGA input 
TCGA_input = loadRData(paste(PATH,'/LIRICS/organized data/TCGA/TCGA_LIRICS.Rdata',sep=''))
TCGA_surv_input =  loadRData(paste(PATH,'/LIRICS/organized data/TCGA/TCGA_LIRICS_SURV.Rdata',sep=''))

## load LIRICS relevant interactions 
lirics_relevant = readRDS(paste(PATH,'/LIRICS/LIRICS database/lirics_skcm_database_interactions.RDS',sep=''))

## SOCIAL ## 

## load SOCIAL output
ICB_sc_input = readRDS(paste(PATH,'/SOCIAL/output/Jerby-Arnon et al. Cell 2018/livnat2018_TPM_400pseudopat_40ds_SOCIAL.rds',sep=''))

## IRIS ## 

## load selected features from IRIS
feature_RDI = loadRData(paste(PATH,'/IRIS/results/IRIS_RDI_output.Rdata',sep=''))
feature_RUI = loadRData(paste(PATH,'/IRIS/results/IRIS_RUI_output.Rdata',sep=''))

feature_RDI_mc = loadRData(paste(PATH,'/IRIS/results/IRIS_RDI_output_mc.Rdata',sep='')) #testing features for monotherapy combination therapy seperate
feature_RUI_mc = loadRData(paste(PATH,'/IRIS/results/IRIS_RUI_output_mc.Rdata',sep='')) #testing features for monotherapy combination therapy seperate


#load SPECIAL output 
thrane_output = readRDS(paste(PATH,'/3. SOCIAL & SPECIAL/SPECIAL/thrane 2018 st/output/thrane2018_SPECIAL_output.RDS',sep=''))
biermann_output = readRDS(paste(PATH,'/3. SOCIAL & SPECIAL/SPECIAL/biermann 2022 st/output/biermann22_SPECIAL_output.RDS',sep=''))

## Benchmark ##

#TIDE
TIDE_AUC = readRDS(paste(PATH,'/benchmark/TIDE/tide_auc.RDS',sep=''))
TIDE_AUC_supp= readRDS(paste(PATH,'/benchmark/TIDE/tide_auc_supp.RDS',sep='')) #mono comb therapy specific 
TIDE_score = readRDS(paste(PATH,'/benchmark/TIDE/tide_input.RDS',sep=''))

#signatures except TIDE (including IMPRES)
sig_AUC = readRDS(paste(PATH,'/benchmark/signatures/sig_auc.RDS',sep=''))
sig_score = readRDS(paste(PATH,'/benchmark/signatures/sig_input.RDS',sep=''))
sig_thrane = readRDS('/benchmark/signatures/sig_input_spatial.RDS') #CHANGE Directory

#IMPRES only AUC
IMPRES_AUC_supp= readRDS(paste(PATH,'/benchmark/signatures/impres_auc_supp.RDS',sep='')) #mono comb therapy specific 

## Miscellaneous ##
pcg = fread(paste(PATH,'/miscellaneous/coding_genes.txt',sep=''), header=F) #protein coding genes

