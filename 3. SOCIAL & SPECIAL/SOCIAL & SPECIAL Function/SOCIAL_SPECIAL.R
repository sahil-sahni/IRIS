## SOCIAL & SPECIAL
library(tidyverse)
library(Matrix)
library(parallel)
## ------- Functions for SOCIAL -------

## Step 1 Functions:

#SOCIAL.query_LRdb(): A function to extract ligand receptor pairs that are found within our sc data. 
#db: a data frame of ligand receptor gene pairs, column V1 is the ligand, column V2 is the receptor (complex with up-to 3 genes).
#expr: expression matrix (log normalized counts), rows are genes, columns are cells.  

#output:
#pairs: a data frame of ligand receptor gene pairs, column V1 is the ligand, column V2 is the receptor
#it is of the following format:
#V1 V2
#"A2M"    "LRP1"  
#"ACE"    "BDKRB2"
#"ADAM10" "AXL"   
#"ADAM10" "EPHA3" 
#"ADAM12" "ITGB1" 
#"ADAM12" "SDC4" 

SOCIAL.query_LRdb <- function(db, expr){
  
  sc_gene = toupper(rownames(expr))
  
  pairs = toupper(as.matrix(db))
  pairs = unique(pairs)
  
  pairs_check = pairs %>% as.data.frame(.) %>% mutate(receptor=str_split_fixed(pairs[,2], ";", 3)) %>% mutate_all(na_if, "") %>% as.matrix(.) 
  
  ind1 = (pairs_check[,1] %in% sc_gene); 
  ind2 = (pairs_check[,2] %in% sc_gene);
  ind3 = (pairs_check[,3] %in% sc_gene)|(pairs_check[,3] %in% c(""));
  ind4 = (pairs_check[,4] %in% sc_gene)|(pairs_check[,4] %in% c(""));
  
  pairs = pairs[ind1&ind2&ind3&ind4,]
  
  return(pairs)
}

## Step 2 & 3 Functions:

#SOCIAL.calculate_interaction_score(): A function to infer cell-2-cell interactions from sc data.
#expr: expression matrix (log normalized counts), rows are genes, columns are cells.  
#Set rownames to gene symbols and column names to cell ids


#ct_map: a dataframe describing the mapping of each cell to individual and cell type. it is of the following format:
#Cell cell_type samples
#         cy79_p4_CD45_neg_PDL1_neg_E11_S1115_comb       EVT   Mel79
#  cy79_p1_CD45_neg_PDL1_pos_AS_C1_R1_F07_S67_comb       EVT   Mel79
#  cy79_p1_CD45_neg_PDL1_pos_AS_C1_R1_G01_S73_comb       EVT   Mel79
# cy79_p1_CD45_neg_PDL1_neg_AS_C4_R1_D09_S141_comb       EVT   Mel79
# cy79_p1_CD45_neg_PDL1_neg_AS_C4_R1_D06_S138_comb       EVT   Mel79
#                    cy53_1_CD45_neg_C06_S318_comb       EVT   Mel53

#pairs: a data frame of ligand receptor gene pairs, column V1 is the ligand, column V2 is the receptor
#it is of the following format:
#V1 V2
#"A2M"    "LRP1"  
#"ACE"    "BDKRB2"
#"ADAM10" "AXL"   
#"ADAM10" "EPHA3" 
#"ADAM12" "ITGB1" 
#"ADAM12" "SDC4"  

#n_iterations: number of bootstrap iterations to infer the null distribution

SOCIAL.calculate_interaction_score <- function(expr,ct_map, pairs, n_iterations = 100){ #for spatial data
  
  LR_interest = c(pairs$ligand, pairs$receptor) %>% strsplit(., ";") %>% unlist(.) %>% unique(.) # LR of interest to speed of data aggregation
  
  ctypes = unique(ct_map$cell_type)
  ctype_groupings = t(combn(ctypes,2))
  ctype_groupings = rbind(ctype_groupings,ctype_groupings[,c(2,1)])
  
  #compute real scores (Step 2)
  summary_expr = aggregate.data.frame(as.data.frame(t(expr[LR_interest,])), by = list(ct_map$cell_type), FUN = mean, na.rm = T)
  rownames(summary_expr) = summary_expr$Group.1
  summary_expr = summary_expr[,-1]
  summary_expr = t(summary_expr)
  
  summary_expr2 = aggregate.data.frame(as.data.frame(t(expr[LR_interest,]) > 0), by = list(ct_map$cell_type), FUN = mean, na.rm = T)
  rownames(summary_expr2) = summary_expr2$Group.1
  summary_expr2 = summary_expr2[,-1]
  summary_expr2 = t(summary_expr2)
  
  real_scores = sapply(1:nrow(ctype_groupings), function(j) { # compute interaction scores
    ct1 = ctype_groupings[j,1]
    ct2 = ctype_groupings[j,2]
    scores = suppressMessages(t(apply(pairs,1,function(y){
      g1 = strsplit(y[1],split = ";")[[1]]
      g2 = strsplit(y[2],split = ";")[[1]]
      score = c()
      ligand_exp = NULL
      receptor_exp = NULL
      if(length(g1) > 1)
        ligand_exp = min(as.numeric(summary_expr[g1,ct1])) 
      else
        ligand_exp = as.numeric(summary_expr[g1,ct1])
      if(length(g2) > 1)
        receptor_exp = min(as.numeric(summary_expr[g2,ct1]))
      else
        receptor_exp = as.numeric(summary_expr[g2,ct1])
      
      #interaction score
      sc = ligand_exp*receptor_exp 
      return(sc)
    })))
    return(scores)
  })
  
  ligand_scores = sapply(1:nrow(ctype_groupings), function(j) { # compute mean ligand expression
    ct1 = ctype_groupings[j,1]
    ct2 = ctype_groupings[j,2]
    scores = suppressMessages(t(apply(pairs,1,function(y){
      g1 = strsplit(y[1],split = ";")[[1]]
      score = c()
      ligand_exp = NULL
      if(length(g1) > 1)
        ligand_exp = min(as.numeric(summary_expr[g1,ct1])) 
      else
        ligand_exp = as.numeric(summary_expr[g1,ct1])
      
      #ligand expression
      sc = ligand_exp
      return(sc)
    })))
    return(scores)
  })
  
  receptor_scores = sapply(1:nrow(ctype_groupings), function(j) { # compute mean receptor expression
    ct1 = ctype_groupings[j,1]
    ct2 = ctype_groupings[j,2]
    scores = suppressMessages(t(apply(pairs,1,function(y){
      g2 = strsplit(y[2],split = ";")[[1]]
      score = c()
      receptor_exp = NULL
      if(length(g2) > 1)
        receptor_exp = min(as.numeric(summary_expr[g2,ct1]))
      else
        receptor_exp = as.numeric(summary_expr[g2,ct1])
      
      #receptor expression
      sc = receptor_exp 
      return(sc)
    })))
    return(scores)
  })
  
  rownames(real_scores) = sapply(1:nrow(pairs), function(i) paste(pairs[i,1],"_",pairs[i,2],sep = ""))
  colnames(real_scores) = sapply(1:nrow(ctype_groupings), function(i) paste(ctype_groupings[i,1],"_",ctype_groupings[i,2],sep = ""))
  
  rownames(ligand_scores) = sapply(1:nrow(pairs), function(i) paste(pairs[i,1],"_",pairs[i,2],sep = ""))
  colnames(ligand_scores) = sapply(1:nrow(ctype_groupings), function(i) paste(ctype_groupings[i,1],"_",ctype_groupings[i,2],sep = ""))
  
  rownames(receptor_scores) = sapply(1:nrow(pairs), function(i) paste(pairs[i,1],"_",pairs[i,2],sep = ""))
  colnames(receptor_scores) = sapply(1:nrow(ctype_groupings), function(i) paste(ctype_groupings[i,1],"_",ctype_groupings[i,2],sep = ""))
  
  #compute null scores (Step 3)
  null_scores = lapply(1:n_iterations, function(i) {
    start = Sys.time()
    usam = unique(ct_map$samples)
    null_labels = ct_map$cell_type
    for(samp in usam){ 
      null_labels[ct_map$samples == samp] = sample(ct_map$cell_type[ct_map$samples == samp]) #shuffle cell labels
    }
    null_map = data.frame(cell = ct_map$Cell, cell_type = sample(null_labels))
    summary_expr = aggregate.data.frame(as.data.frame(t(expr[LR_interest,])), by = list(null_map$cell_type), FUN = mean, na.rm = T)
    rownames(summary_expr) = summary_expr$Group.1
    summary_expr = summary_expr[,-1]
    summary_expr = t(summary_expr)
    
    summary_expr2 = aggregate.data.frame(as.data.frame(t(expr[LR_interest,]) > 0), by = list(null_map$cell_type), FUN = mean, na.rm = T)
    rownames(summary_expr2) = summary_expr2$Group.1
    summary_expr2 = summary_expr2[,-1]
    summary_expr2 = t(summary_expr2)
    null_score = sapply(1:nrow(ctype_groupings), function(j) {
      ct1 = ctype_groupings[j,1]
      ct2 = ctype_groupings[j,2]
      scores = suppressMessages(t(apply(pairs,1,function(y){
        g1 = strsplit(y[1],split = ";")[[1]]
        g2 = strsplit(y[2],split = ";")[[1]]
        score = c()
        ligand_exp = NULL
        receptor_exp = NULL
        if(length(g1) > 1)
          ligand_exp = min(as.numeric(summary_expr[g1,ct1])) 
        else
          ligand_exp = as.numeric(summary_expr[g1,ct1])
        if(length(g2) > 1)
          receptor_exp = min(as.numeric(summary_expr[g2,ct1]))
        else
          receptor_exp = as.numeric(summary_expr[g2,ct1]) 
        
        #null score
        sc = ligand_exp*receptor_exp 
        return(sc)
      })))
      return(scores)
    })
    print(Sys.time() - start)
    return(null_score)
  })
  
  #compute empirical p-value (Step 3)
  pvals = matrix(0,nrow = nrow(real_scores), ncol = ncol(real_scores))
  for(i in 1:n_iterations){
    pvals = pvals + (null_scores[[i]] >= real_scores) #assumption is >= so that 0 >= 0 if considered true
  }
  pvals = pvals/n_iterations
  rownames(pvals) = sapply(1:nrow(pairs), function(i) paste(pairs[i,1],"_",pairs[i,2],sep = ""))
  colnames(pvals) = sapply(1:nrow(ctype_groupings), function(i) paste(ctype_groupings[i,1],"_",ctype_groupings[i,2],sep = ""))
  
  return(list(interaction_scores = real_scores,p.value = pvals, ligand_exp=ligand_scores, receptor_exp=receptor_scores))
}

## Step 4 Functions:

#SOCIAL.infer_activity(): A function to infer ligand receptor interaction activity from SOCIAL object.
#input: SOCIAL object outputted from SOCIAL.cis_rslurm() or SOCIAL.calculate_interaction_score.
#p: empirical p-value cutoff (what empirical p value is considered significantly active)
#median: whether activity should implement ligand and receptor median expression (above median = active)

SOCIAL.infer_activity <- function(input, p_cutoff=0.05, median_cutoff=T){
  
  # empirical p value (active if p < 0.05)
  input_p = input$empiricalp %>% set_rownames(.$int) %>% .[,-1] %>% as.matrix(.) 
  input_p = ifelse(input_p < p_cutoff, 1, 0) %>% ifelse(is.na(.), 0, .)
  
  if(median_cutoff==T){
    # ligand median expression (active if mean ligand expr > median)
    input_l = input$ligandexp %>% set_rownames(.$int) %>% .[,-1] %>% as.matrix(.) 
    input_l = apply(input_l, 1, function(x){ifelse(x > median(x, na.rm=T), 1,0) }) %>% t(.) %>% ifelse(is.na(.), 0, .)
    
    # receptor median expression (active if mean receptor expr > median)
    input_r = input$receptorexp %>% set_rownames(.$int) %>% .[,-1] %>% as.matrix(.) 
    input_r = apply(input_r, 1, function(x){ifelse(x > median(x, na.rm=T), 1,0) }) %>% t(.) %>% ifelse(is.na(.), 0, .)
    
    # create a 3D array
    exp_array <- abind(input_p, input_l, input_r, along=3)
    dimnames(exp_array)[[3]] = list('cci','ligand','receptor')
    dim(exp_array)
    
    # sum across 3 array cci, ligand, receptor (1+1+1 = 3 if active)
    sum_array <- apply(exp_array, c(1, 2), sum)
    
    # spot score where sum == 3 is 1 (active); others is 0 (inactive)
    output_3 = ifelse(sum_array == 3, 1, 0) %>% t(.)
    return(output_3)
  }
  
  output_1 = input_p
  return(output_1)
}

## Additional Functions:

#SOCIAL.cis_rslurm(): A function to infer cell-2-cell interactions (SOCIAL steps 2-3) from sc data (for each sample) using rslurm (ideal for when multiple samples are available).

#expr: expression matrix (log normalized counts), rows are genes, columns are cells.  
#Set rownames to gene symbols and column names to cell ids

#ct_map: a dataframe describing the mapping of each cell to individual and cell type. it is of the following format:
#Cell cell_type samples
#         cy79_p4_CD45_neg_PDL1_neg_E11_S1115_comb       EVT   Mel79
#  cy79_p1_CD45_neg_PDL1_pos_AS_C1_R1_F07_S67_comb       EVT   Mel79
#  cy79_p1_CD45_neg_PDL1_pos_AS_C1_R1_G01_S73_comb       EVT   Mel79
# cy79_p1_CD45_neg_PDL1_neg_AS_C4_R1_D09_S141_comb       EVT   Mel79
# cy79_p1_CD45_neg_PDL1_neg_AS_C4_R1_D06_S138_comb       EVT   Mel79
#                    cy53_1_CD45_neg_C06_S318_comb       EVT   Mel53

#pairs: a data frame of ligand receptor gene pairs, column V1 is the ligand, column V2 is the receptor
#it is of the following format:
#V1 V2
#"A2M"    "LRP1"  
#"ACE"    "BDKRB2"
#"ADAM10" "AXL"   
#"ADAM10" "EPHA3" 
#"ADAM12" "ITGB1" 
#"ADAM12" "SDC4"  

#n_iterations: number of bootstrap iterations to infer the null distribution

#path_social: path to SOCIAL function
#path_sc: path to sc expr matrix

SOCIAL.cis_rslurm <- function(expr,ct_map, pairs, n_iterations = 100, path_social){
  
  ## load SOCIAL functions
  source(path_social)
  
  ## create jobname 
  jobname = paste('SOCIAL', Sys.time())
  print(jobname)
  
  ## create directory to save intermediate input
  dir = getwd()
  dir.create(paste(dir,'/',jobname,'/',sep=""))
  temp_path_cis = paste(dir,'/',jobname,'/',sep="")
  
  ## save sc expr files as sparse matrix
  sparse_expr=Matrix(expr,sparse = T)
  sparse_path = paste(temp_path_cis,'sparse_expr.RDS',sep='')
  saveRDS(sparse_expr, file=sparse_path)
  
  ## generate and save intermediate ct_map files
  ct_map_input = generate_ctmap_per_sample(ct_map, expr)
  ct_map_path = paste(temp_path_cis,'intermediate_ctmap.RDS',sep='')
  saveRDS(ct_map_input, file=ct_map_path)
  
  ## save pairs files
  pairs_path = paste(temp_path_cis,'LR_pairs_prefiltered.csv',sep='')
  write.csv(pairs, pairs_path) 
  
  ## generate rslurm data.frame
  run = 1:length(unique(ct_map$samples))
  samples = unique(ct_map$samples)
  expr_path = sparse_path
  ct_map_path = ct_map_path
  pairs_path = pairs_path
  path_social = path_social
  n_iterations= n_iterations
  
  rslurm_df = data.frame(run=run, samples=samples, expr_path=expr_path, ct_map_path=ct_map_path, pairs_path=pairs_path, path_social=path_social, n_iterations=n_iterations)

  ## run rslurm
  sjob = slurm_apply(run_cis_rslurm, rslurm_df, jobname=jobname, nodes=400, cpus_per_node=1,submit=TRUE, slurm_options=list(time='24:00:00', mem='20g')) #CHANGE HYPERPARAMETER AS NEEDED
  profile = get_slurm_out(sjob, outtype = 'raw', wait = TRUE) #SOCIAL results
  names(profile) = rslurm_df$samples
  
  ## aggregate results
  pval_list = lapply(names(profile), function(x){ 
    m_pv = melt(profile[[x]]$p.value)
    colnames(m_pv) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_pv = m_pv %>% mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)] 
    return(m_pv)
  })
  
  score_list = lapply(names(profile), function(x){
    m_sc = melt(profile[[x]]$interaction_scores)
    colnames(m_sc) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_sc = m_sc %>% mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)]
    return(m_sc)
  })
  
  ligand_list = lapply(names(profile), function(x){
    m_pv = melt(profile[[x]]$ligand_exp)
    colnames(m_pv) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_pv = m_pv %>% mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)] 
    return(m_pv)
  })
  
  receptor_list = lapply(names(profile), function(x){
    m_sc = melt(profile[[x]]$receptor_exp)
    colnames(m_sc) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_sc = m_sc %>% mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)]
    return(m_sc)
  })
  
  pval_df = pval_list %>% reduce(full_join, by='int')
  scr_df = score_list %>% reduce(full_join, by='int')
  ligand_df = ligand_list %>% reduce(full_join, by='int')
  receptor_df = receptor_list %>% reduce(full_join, by='int')
  
  unlink(temp_path_cis, recursive = TRUE) #remove temporary path
  
  social_profile=list(empiricalp=pval_df, interactionscore=scr_df, ligandexp=ligand_df, receptorexp=receptor_df)
  return(social_profile)
}

## Internal Functions:

#generate_ctmap_per_sample(): A internal function to generates ct_map for each sample to load faster for generating interaction score using cis_rslurm().

#ct_map: a dataframe describing the mapping of each cell to individual and cell type. it is of the following format:
#Cell cell_type samples
#         cy79_p4_CD45_neg_PDL1_neg_E11_S1115_comb       EVT   Mel79
#  cy79_p1_CD45_neg_PDL1_pos_AS_C1_R1_F07_S67_comb       EVT   Mel79
#  cy79_p1_CD45_neg_PDL1_pos_AS_C1_R1_G01_S73_comb       EVT   Mel79
# cy79_p1_CD45_neg_PDL1_neg_AS_C4_R1_D09_S141_comb       EVT   Mel79
# cy79_p1_CD45_neg_PDL1_neg_AS_C4_R1_D06_S138_comb       EVT   Mel79
#                    cy53_1_CD45_neg_C06_S318_comb       EVT   Mel53

#expr: expression matrix (log normalized counts), rows are genes, columns are cells. columns is non-subsetted by samples!

generate_ctmap_per_sample <- function(ct_map, expr){
  samples = unique(ct_map$samples)
  
  ct_map_list = mclapply(1:length(samples), function(i){ #generates ct_map files for each sample
    sample_i = samples[i]
    ct_map_sample = subset(ct_map, (Cell %in% colnames(expr)) & (samples == sample_i)) #ct_map
    return(ct_map_sample)
  })
  
  names(ct_map_list) = samples
  return(ct_map_list)
}

#run_cis_rslurm(): A internal function to execute calculate interaction score on rslurm efficiently

run_cis_rslurm <- function(run, samples, expr_path, ct_map_path, pairs_path, path_social, n_iterations){
  
  #load functions
  source(path_social)
  
  #CODE START
  set.seed(run+1)
  
  #load sc data
  expr = readRDS(expr_path) %>% as.matrix(.)
  print('expr loaded')
  
  #load filtered LR pairs
  pairs = read.csv(pairs_path) %>% dplyr::select(ligand,receptor)
  print('pairs loaded')
  
  #load ct_map
  load_ctmap = readRDS(ct_map_path)
  ct_map = load_ctmap[[samples]] %>% dplyr::select(Cell, cell_type, samples)
  print('ctmap loaded')
  
  #filter expr (assuming duplicate cell IDs):
  ct_map = ct_map %>% mutate(uID = paste(Cell,1:nrow(ct_map), sep="nci")) #append character to help prevent duplicates
  id_df = data.frame(Cell=ct_map$Cell, uID=ct_map$uID) 
  expr_filter = merge(id_df, t(expr), by.x=1, by.y=0) %>% magrittr::set_rownames(.$uID) %>% .[,-1:-2] %>% as.matrix(.) %>% t(.)
  ct_map = ct_map %>% mutate(Cell = uID) %>% select(-uID)
  expr_filter = expr_filter[,ct_map$Cell] #reorder
  colnames(expr_filter)
  
  SOCIAL_res = SOCIAL.calculate_interaction_score(expr_filter,ct_map, pairs, n_iterations)
  print('SOCIAL Step 2&3 Complete')
  
  return(SOCIAL_res)
}

## ------- Functions for SPECIAL -------

## Main Function:
#SPECIAL(): function to execute steps 1-2 of SPECIAL
#
#expr: expression matrix (log normalized counts), rows are genes, columns are cells.  
#Set rownames to gene symbols and column names to cell ids

#loc: a dataframe describing the mapping of each cell to spatial coordinate and cell type (output of CytoSPACE. it is of the following format:
#   UniqueCID: unique sc barcode (generated by CytoSPACE)
#   OriginalCID:  original sc barcode from study
#   CellType: cell type of individual single-cell
#   SpotID: barcode or coordinate of the spatial location of the aligned single-cell
#   row: y coordinate of aligned single-cell
#   col: x coordinate of aligned single-cell

#pairs: a data frame of ligand receptor gene pairs, column V1 is the ligand, column V2 is the receptor
#it is of the following format:
#V1 V2
#"A2M"    "LRP1"  
#"ACE"    "BDKRB2"
#"ADAM10" "AXL"   
#"ADAM10" "EPHA3" 
#"ADAM12" "ITGB1" 
#"ADAM12" "SDC4"  

#n_iterations: number of bootstrap iterations to infer the null distribution
#distance: distance of spatial region (for bulk it is the radius; for SlideSeq it is the diameter)
#platform: 'slideseqv2' or 'bulk'
#puck_diameter: diameter (in unit "x") of the SlideSeqV2 puck (for SlideSeqV2 ONLY) (optional argument)

#path_social: path to SOCIAL function
SPECIAL <- function(expr, pairs, loc, name, n_iterations, path_special, distance, platform, puck_diameter=NA){
  print('running SPECIAL')
  
  ## source SPECIAL
  source(path_special)
  
  ## initialize seed
  set.seed(20892)
  
  ## create directory to save intermediate input
  DIR = getwd()
  dir.create(paste(DIR,'/temporary/',sep=""))
  path_1 = paste(DIR,'/temporary/',sep="")
  dir.create(paste(path_1,name,'/',sep="")) 
  temp_path_cis = paste(path_1, name, '/', sep='') #path to temporary folder
  
  ## run step 1 of SPECIAL
  
  # if legacy or visium
  if (toupper(platform) == 'BULK'){
    print('input bulk')
    
    region = SPECIAL.sliding_window(loc, distance)
  }
  
  # if slideseqv2 
  if (toupper(platform) == 'SLIDESEQV2'){
    print('input SLIDESEQV2')
    
    region = SPECIAL.cluster_region(loc, distance, puck_diameter)
  }
  
  ## calculate cell fraction
  
  cellcount = region %>% 
    group_by(`samples`, cell_type) %>%
    dplyr::summarise(count=n()) %>% 
    pivot_wider(., names_from=cell_type, values_from=count) %>% as.data.frame(.) %>%
    set_rownames(.$samples) %>% .[,-1] %>% replace(., is.na(.), 0) 
  cellfraction = cellcount/rowSums(cellcount)  
  
  ## remove regions ("samples" ~ to seamlessly apply to SOCIAL) and cellfraction with only one cell type 
  removed_regions <- names(which(apply(cellfraction == 1, 1, any)))
  
  cellfraction = cellfraction[!(rownames(cellfraction) %in% removed_regions),] 
  region = subset(region, !(samples %in% removed_regions))
  
  ## generate intermediate ct_map per region ("samples")
  ct_map_input = generate_ctmap_per_sample(region, expr)
  ct_map_path = paste(temp_path_cis,'intermediate_ctmap.RDS',sep='')
  saveRDS(ct_map_input, file=ct_map_path)
  
  ## save sc expr files as sparse matrix
  sparse_expr=Matrix(expr,sparse = T)
  sparse_path = paste(temp_path_cis,'sparse_expr.RDS',sep='')
  saveRDS(sparse_expr, file=sparse_path)
  
  ## save pairs files
  pairs_path = paste(temp_path_cis,'LR_pairs_prefiltered.csv',sep='')
  write.csv(pairs, pairs_path) 
  
  ## generate rslurm data.frame
  run = 1:length(unique(region$samples))
  samples = unique(region$samples)
  expr_path = sparse_path
  ct_map_path = ct_map_path
  pairs_path = pairs_path
  path_social = path_special
  n_iterations = n_iterations
  
  rslurm_df = data.frame(run=run, samples=samples, expr_path=expr_path, ct_map_path=ct_map_path, pairs_path=pairs_path, path_social=path_social, n_iterations=n_iterations)
  
  ## run rslurm (SOCIAL Step 2-3)
  sjob = slurm_apply(run_cis_rslurm, rslurm_df, jobname=name, nodes=nrow(rslurm_df), cpus_per_node=8,submit=TRUE, slurm_options=list(time='72:00:00', mem='280g', partition='norm,ccr'), preschedule_cores=FALSE) 
  profile = get_slurm_out(sjob, outtype = 'raw', wait = TRUE)
  cleanup_files(sjob) #remove slurm folder
  names(profile) = rslurm_df$samples
  
  ## creates data.frame of all possible ligand_receptor pairs x samples; values are empirical p
  pval_list = lapply(names(profile), function(x){
    m_pv = reshape2::melt(profile[[x]]$p.value)
    colnames(m_pv) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_pv = m_pv %>% mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)] 
    return(m_pv)
  })
  
  pval_df = pval_list %>% reduce(full_join, by='int')
  
  ## creates data.frame of all possible ligand_receptor pairs x samples; values are ligand exp
  ligand_list = lapply(names(profile), function(x){
    m_pv = reshape2::melt(profile[[x]]$ligand_exp)
    colnames(m_pv) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_pv = m_pv %>% mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)] 
    return(m_pv)
  })
  
  ligand_df = ligand_list %>% reduce(full_join, by='int')
  
  ## creates data.frame of all possible ligand_receptor pairs x samples; values are receptor exp
  receptor_list = lapply(names(profile), function(x){
    m_pv = reshape2::melt(profile[[x]]$receptor_exp)
    colnames(m_pv) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_pv = m_pv %>% mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)] 
    return(m_pv)
  })
  
  receptor_df = receptor_list %>% reduce(full_join, by='int')
  
  ## delete created directory
  #unlink(temp_path_cis, recursive = TRUE)
  
  ## return objects
  return(list(name = name, cellfraction=cellfraction, empiricalp=pval_df, removed_regions=removed_regions, ligandexp=ligand_df, receptorexp=receptor_df, cellcount=cellcount))
}


## Step 1 Functions:

#SPECIAL.cluster_region(): function to partiion spatial transcriptomics slide to physiological relevant ligand-receptor interaction regions
#loc: a dataframe describing the mapping of each cell to spatial coordinate and cell type (output of CytoSPACE. it is of the following format:
#   UniqueCID: unique sc barcode (generated by CytoSPACE)
#   OriginalCID:  original sc barcode from study
#   CellType: cell type of individual single-cell
#   SpotID: barcode or coordinate of the spatial location of the aligned single-cell
#   row: y coordinate of aligned single-cell
#   col: x coordinate of aligned single-cell
#cluster_diameter: diameter (in unit "x") of the individual regions
#puck_diameter: diameter (in unit "x") of the SlideSeqV2 puck

SPECIAL.cluster_region <- function(loc, cluster_diameter, puck_diameter){
  
  # create puck coordinate map
  map = as.data.frame(loc) %>% mutate(x=row, y=col) 
  map = map %>% mutate(x = (x-min(x))/max(x)*puck_diameter, y = (y-min(y))/max(y)*puck_diameter) #rescale to puck dimension
  coord = as.matrix(map[,c('x','y')]) #create matrix
  
  # initialize 
  puck_radius = puck_diameter/2
  cluster_radius = cluster_diameter/2
  centers = as.integer((pi*puck_radius^2)/(pi*cluster_radius^2)) #number of spatial spots that can fit in the puck
  
  # calculate kmean clustering with ideal centers
  set.seed(20892)
  kmeans_result = kmeans(coord, centers=centers)
  map$cluster = kmeans_result$cluster
  
  # cluster center data.frame
  cluster_coord = data.frame(cluster = 1:centers, kmeans_result$centers) %>% mutate(regions=paste(round(x,4),round(y,4),sep='x'))
  
  # merge data.frame and filter
  res = merge(map, cluster_coord, by.x='cluster', by.y='cluster') %>% select(Cell=OriginalCID, cell_type=CellType, samples=regions)
  return(res)
}

#SPECIAL.sliding_window(): function to partiion spatial transcriptomics slide to physiological relevant ligand-receptor interaction regions
#loc: a dataframe describing the mapping of each cell to spatial coordinate and cell type (output of CytoSPACE. it is of the following format:
#   UniqueCID: unique sc barcode (generated by CytoSPACE)
#   OriginalCID:  original sc barcode from study
#   CellType: cell type of individual single-cell
#   SpotID: barcode or coordinate of the spatial location of the aligned single-cell
#   row: y coordinate of aligned single-cell
#   col: x coordinate of aligned single-cell
#radius: radius (in unit "x") of the individual region

SPECIAL.sliding_window <- function(loc, radius=1){
  map = data.frame(x=loc$row, y=loc$col) %>% distinct(.)
  
  windows = lapply(1:nrow(map), function(o){
    origin = map[o,]
    
    window = subset(loc, sqrt((row - origin$x)^2 + (col - origin$y)^2) <= radius) #euclidean distance
    window$window = paste(origin$x,'x',origin$y,sep='')
    
    window = window %>% select(Cell=OriginalCID, cell_type=CellType, samples=window)
    return(window)
  })
  windows = do.call(rbind,windows)
  return(windows)
}

## Additional Functions:

#SPECIAL.cohort(): A function to infer cell-2-cell interactions (SPECIAL steps 1-2) from spatial data (for each spatial sample) using rslurm (ideal for when multiple samples are available).
#
#path_special: path to SPECIAL function
#
#path_input: path to SPECIAL organized input; two layer list with the outer layer corresponding to each spatial sample slide; while inner layer contains the following objects:
#   name: name of spatial slide sample
#   cci: a data frame of ligand receptor gene pairs, column V1 is the ligand, column V2 is the receptor; this step is done by running (SOCIAL.query_LRdb) on expr data
#        it is of the following format:
#        V1 V2
#        "A2M"    "LRP1"  
#        "ACE"    "BDKRB2"
#        "ADAM10" "AXL"   
#        "ADAM10" "EPHA3" 
#        "ADAM12" "ITGB1" 
#        "ADAM12" "SDC4" 
#   loc: data.frame corresponding to the individual single-cell locations within an spatial transcriptomics slide (CytoSPACE output)
#   contains the following columns:
#     UniqueCID: unique sc barcode (generated by CytoSPACE)
#     OriginalCID:  original sc barcode from study
#     CellType: cell type of individual single-cell
#     SpotID: barcode or coordinate of the spatial location of the aligned single-cell
#     row: y coordinate of aligned single-cell
#     col: x coordinate of aligned single-cell
#   exp: expression matrix (log normalized counts), rows are genes, columns are cells.  
#   Set rownames to gene symbols and column names to original cell ids
#
#n_iterations: number of bootstrap iterations to infer the null distribution (for SOCIAL)
#distance: distance of spatial region (for bulk it is the radius; for SlideSeq it is the diameter)
#platform: 'SlideSeqv2' or 'Bulk'
#puck_diameter: for SlideSeq

SPECIAL.cohort <- function(path_special, path_input, n_iterations, distance, platform, puck_diameter=NA){
  
  print('running SPECIAL.cohort')
  
  #load SPECIAL and its input
  source(path_special)
  special_input = readRDS(path_input)
  
  #create data.frame to run rslurm across different slides/samples
  rslurm_df = data.frame(index=1:length(special_input), path_input=path_input, path_special = path_special, n_iterations=n_iterations, distance=distance, platform=platform, puck_diameter=puck_diameter)
  jobname = paste('per','index',Sys.time(), sep='_')
  
  #run rslurm 
  sjob_out = slurm_apply(SPECIAL.per_sample, rslurm_df, jobname=jobname, nodes=nrow(rslurm_df), cpus_per_node=8,submit=TRUE, slurm_options=list(time='72:00:00', mem='10g', partition='norm,ccr'), preschedule_cores=FALSE) 
  output = get_slurm_out(sjob_out, outtype = 'raw', wait = TRUE)
  #cleanup_files(sjob_out)
  
  return(output)
}

## Internal Functions:

#SPECIAL.per_sample(): A internal function to execute SPECIAL on each sample. Can be expanded in the future to add Null distribution.

SPECIAL.per_sample <- function(index, path_input, path_special, n_iterations, distance, platform, puck_diameter=NA){
  print('running SPECIAL.per_sample')
  
  #load SPECIAL
  source(path_special)

  #data.frame for rslurm
  rslurm_df = data.frame(index=index, path_input=path_input, path_special=path_special, n_iterations=n_iterations, distance=distance, platform=platform, puck_diameter=puck_diameter)
  
  #run rslurm per slide
  sjob_in = slurm_apply(run_special_rslurm, rslurm_df, jobname=paste('index',index,sep="_"), nodes=nrow(rslurm_df), cpus_per_node=10,submit=TRUE, slurm_options=list(time='72:00:00', mem='40g', partition='norm,ccr'), preschedule_cores=FALSE) #2 before 
  output = get_slurm_out(sjob_in, outtype = 'raw', wait = TRUE)
  #cleanup_files(sjob_in)
  
  #outputs
  return(output)
}

#run_special_rslurm(): A internal function to execute SPECIAL on rslurm

run_special_rslurm <- function(index, path_input, path_special, name, n_iterations, distance, platform, puck_diameter=NA){
  print('running run_special_rslurm')
  
  #load special functions
  source(path_special)
  
  #load input path
  input = readRDS(path_input)
  
  #load indexed sample
  indexed_sample = input[[index]]
  expr = indexed_sample$expr
  loc = indexed_sample$loc
  name = indexed_sample$name
  pairs = indexed_sample$pairs
  distance = distance
  platform = platform
  puck_diameter = puck_diameter
  
  #run SPECIAL step 1 and 2
  results = SPECIAL(expr, pairs, loc, name, n_iterations, path_special, distance, platform, puck_diameter)
  return(results)
}