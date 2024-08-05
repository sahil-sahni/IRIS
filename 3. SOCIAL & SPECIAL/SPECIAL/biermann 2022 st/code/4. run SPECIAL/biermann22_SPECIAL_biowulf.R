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

# New file from original submission 07/01/24 SS

## ------- Input -------
setwd("/data/sahnis2/SOCIAL/test_izar") #CHANGE

# path to SPECIAL function (for RSLURM)
path_special = 'SOCIAL_SPECIAL.R' #CHANGE
source(path_special) #CHANGE

# load inputs
path_input = 'biermann2022_SPECIAL_input.RDS'

## ------- Initialize environment -------
## set seed 
set.seed(20892) #NIH zip code

## ------- Run SPECIAL -------
n_iterations = 100
distance = 250
platform = 'slideseqv2'
puck_size = 3000

output = SPECIAL.cohort(path_special, path_input, n_iterations, distance, platform, puck_size)

saveRDS(output, file='biermann22_SPECIAL_output.RDS')
