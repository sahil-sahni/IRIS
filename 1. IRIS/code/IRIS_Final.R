library(tidyverse)
library(magrittr)
library(dplyr)
library(pROC)
library(rBayesianOptimization)
library(parallely)
###

## ------- Major Functions -------
## LOOC_validation(): function used to perform Leave One Out Cohort (LOOC) validation

  ## input:
    ## test_name: string regarding testing cohort name formatted as [cohort]_[pre or on]  (ex. Gide_on, Gide_pre, PUCH_on)
    ## s1_cohort: named list of (data.frame) pre/on-treatment cohort(s) to be used for step1 feature selection; names formatted as [cohort]_[pre/on] (requires at least one pre and one on-rx cohort)
    ## s2_cohort: named list of (data.frame) pre-treatment cohort(s) to be used for step2 feature selection; names formatted as [cohort]_[pre]
    ## direction: string ('RDI' or 'RUI') corresponding to searching for RDIs or RUIs

  ## output:
    ## named list of string of RDIs or RUIs interactions; name correspond to cohort used in step 2 feature selection

LOOC_validation <- function(test_name, s1_cohort, s2_cohort, direction){
  
  print(paste('test cohort:', test_name, sep=" "))
  
  # 1. remove test cohort from step 1 and step 2 feature selection
  test_name = test_name %>% strsplit(., '_') %>% unlist(.) %>% .[[1]] #test cohort name only (excludes pre and on)
    
  s1_name = strsplit(names(s1_cohort),'_') %>% sapply(.,'[[',1) %>% unlist(.) 
  s1_cohort = s1_cohort[which(s1_name != test_name)] #removes testing cohort from step 1 training 
    
  s2_name = strsplit(names(s2_cohort),'_') %>% sapply(.,'[[',1) %>% unlist(.) 
  s2_cohort = s2_cohort[which(s2_name != test_name)] #removes testing cohort from step 2 training
  
  # 2. apply IRIS in a LOOC manner
  IRIS_output =lapply(1:length(s2_cohort), function(i) IRIS(s2_cohort[[i]], names(s2_cohort[i]), s1_cohort, direction)) 
  names(IRIS_output) = names(s2_cohort)
  
  return(IRIS_output)
}

## IRIS(): Immunotherapy Resistant cell-cell Interaction Scanner
  ## input:
    ## s1_cohort: named list of pre/on-treatment cohort(s) to be used for step1 feature selection; names formatted as [cohort]_[pre/on] (requires at least one pre and one on-rx cohort)
    ## s2_cohort: data.frame of pre-treatment  cohort used for step 2 feature selection; names formatted as [cohort]_[pre]
    ## s2_name: name of cohort used for step 2 feature selection names formatted as [cohort]_[pre]
    ## direction: string ('RDI' or 'RUI') corresponding to searching for RDIs or RUIs
  ## output:
    ## vector of selected interactions (Lcell_Rcell_Ligand_Receptor) from IRIS

IRIS <- function(s2_cohort, s2_name, s1_cohort, direction){
  
  # step 1
  
  s1_input = format_s1_input(s2_name, s1_cohort) #format input for step 1
  s1_output = step1(s1_input, direction)
  
  # step 2
  
  fs_mtx = step2(s1_output, s2_cohort[,-1:-2], s2_cohort[,1:2], 500, direction) #matrix of "feature" score from step 2 for each cci from step 1
  feature_score = apply(fs_mtx,2,sum) #sum of all "feature" score for each cci
  
  null_table = data.frame()
  
  for( i in 1:1000){
    set.seed(i)
    random_mtx = t(apply(fs_mtx, 1, sample)) #matrix of shuffled feature score from step 2 for each cci from step 1
    null_score = apply(random_mtx, 2, sum) #sum of all "random" score 
    null_table = rbind(null_table, ifelse(null_score > feature_score, 1, 0))
  }
  
  p_value = (apply(null_table,2,sum))/1000 #calculate empirical p-value for each cci (from step 1) based on feature score generated from step 2 FFS
  sig_index=which(p_value<0.05) #filter for empirical p-values < 0.05
  
  IRIS_output = s1_output[sig_index,]$Interaction #IRIS selected interactions
  
  print('feature selection complete')
  print(length(IRIS_output))
  return(IRIS_output)
}

## ------- Functions for Step 1 -------
## format_s1_input(): function used to ensure training cohort used in step 1 and step 2 feature selection are mutually exclusive and remove testing cohort from training 

  ## input:
    ## s2_name: string regarding testing cohort name formatted as [cohort]_[pre or on]  (ex. Gide_on, Gide_pre, PUCH_on)
    ## s1_input: named list of pre/on-treatment inputs to be used for step1 feature selection; names formatted as [cohort]_[pre/on] (requires at least one pre and one on-rx cohort)

  ## output:
    ## list of 2:
      ## activity_profile: matrix with rows as cci and cols as samples (only pre- and on-NR samples)
      ## meta: data.frame with col: samples; group (1 = ON_NR; 0= PRE)

format_s1_input <- function(s2_name, s1_input){
  
  print(paste('step2 input:',s2_name, sep=' '))
  
  s2_name = s2_name %>% strsplit(., '_') %>% unlist(.) %>% .[[1]] #step 2 cohort name only (i.e. Gide) 
  s1_tp = strsplit(names(s1_input),'_') %>% sapply(.,'[[',2) %>% unlist(.) #step 1 timepoint(s) only (i.e. pre and on)
  
  #1. select all on-treatment cohorts and remove cohort used in step 2
  s1_on = s1_input[which(s1_tp == 'on')] 
  s1_on_name = strsplit(names(s1_on),'_') %>% sapply(.,'[[',1) %>% unlist(.) #step 1 cohort name only
  s1_on = s1_on[which(s1_on_name != s2_name)] #removes cohort used in step 2
  print(paste('step1 on data:',names(s1_on), sep=', '))
  
  #2. merge all on-treatment cohorts and select for samples corresponding to non-responders (0)
  s1_on_nr = do.call(dplyr::bind_rows, s1_on) %>% subset(response==0) #select for non-responder samples 
  row.names(s1_on_nr) = s1_on_nr$sample
  s1_on_nr = s1_on_nr[,-1:-2] %>% t(.) #(rows: cci; col: samples)
  
  #3. select all pre-treatment cohorts and remove cohort used in step 2
  s1_pre = s1_input[which(s1_tp == 'pre')] #select pre-treatment
  s1_pre_name = strsplit(names(s1_pre),'_') %>% sapply(.,'[[',1) %>% unlist(.) #step 1 cohort name only
  s1_pre = s1_pre[which(s1_pre_name != s2_name)] #removes cohort used in step 2
  print(paste('step1 pre data:',names(s1_pre), sep=', '))
  
  #4. merge all pre-treatment cohorts 
  s1_pre_all = do.call(dplyr::bind_rows, s1_pre)
  row.names(s1_pre_all) = s1_pre_all$sample
  s1_pre_all = s1_pre_all[,-1:-2] %>% t(.) #(rows: cci; col: samples)
  
  #5. merge togeather on and pre cohorts as one large activity profile/matrix (rows: cci; col: samples) 
  ap = merge(s1_on_nr, s1_pre_all, by.x=0, by.y=0, all=T)
  row.names(ap) = ap$Row.names
  ap = as.matrix(ap[,-1]) #output 1
  
  #6. create accompanying meta dataframe indicating which samples correspond to on-treatment NR (1) or pre-treatment (0)
  meta = colnames(ap) %>% data.frame(sample=.) %>% mutate(group = ifelse(sample %in% colnames(s1_on_nr),1,0)) #output 2
  
  return(list(activity_profile=ap, meta=meta))
}

## step1(): function used to perform step 1 feature selection (differential activation analysis between pre- and on-NR)

##input:
  ## s1_input: list of 2:
    ## activity_profile: matrix with rows as cci and cols as samples (only pre- and on-NR samples)
    ## meta: data.frame with col including group and sample; group corresponds to either sample belongs to on-rx non-responder (1) or pre-rx (0)
  ## direction: string ('RDI' or 'RUI') corresponding to searching for RDIs or RUIs

## output:
  ## ordered data.frame by ct_FDR (smallest value to largest) with at least three col: 
    ## Interaction (Lcell_Rcell_Ligand_Receptor); 
    ## ct_FDR: FDR per cell pair for each interaction (calculated from p-val regarding differential activation using a fisher.test)
    ## OR: odds-ratio of interaction enriched in NR (>1) or pre-rx (<1) samples

step1 <- function(s1_input, direction){
  
  #1. extract input for fisher test 
  s1_ap = s1_input$activity_profile
  s1_meta = s1_input$meta
  
  #2. run fisher exact test to find differentially activated interactions between pre-treatment and on-treatment non-responder
  fisher_output = mclapply(1:nrow(s1_ap), function(x) step1_fisher(s1_ap[x,], s1_meta)) #fisher test on each interaction individually
  fisher_output = do.call(rbind,fisher_output)
  
  #3. perform multiple hypothesis testing (FDR ~ part 1) for fisher test output 
  fisher_output = data.frame(Interaction=row.names(s1_ap), p=fisher_output$p, OR=fisher_output$OR) %>% 
    mutate(FDR = p.adjust(p,'fdr')) %>% mutate(ct_FDR=NA)
  
  #4. generate/append cell LR interaction data.frame with details of each interaction and their components (Ligand (L), Receptor (R), L cell_type, R cell_type)
  ctLR_df = data.frame(fisher_output$Interaction, str_split_fixed(fisher_output$Interaction, '\\_', 4)) %>%
    set_colnames(.,c("Interaction", "Lcell", "Rcell", "L", "R")) %>% mutate(Lcell_Rcell=paste(Lcell,Rcell,sep='_'))
  
  fisher_output = merge(fisher_output, ctLR_df, by.x=1, by.y=1)
  
  #5. perform multiple hypothesis testing (FDR per celltype pair ~ part 2) for fisher test output 
  cell.pair.counts <- fisher_output$Lcell_Rcell %>% table %>% sort(decreasing = T)       # counts of all cell pairs present
  for (cc in cell.pair.counts %>% names){
    idx <- which(fisher_output$Lcell_Rcell == cc)                                      # indices for the specific cell pair
    fisher_output[idx, "ct_FDR"] <- fisher_output[idx, "p"] %>% 
      p.adjust(method = "fdr")
  }
  
  #6. select for RDI or RUI interactions
  if (direction == 'RDI'){
    print('RDI direction')
    s1_rank = fisher_output %>% subset(ct_FDR < 0.2 & OR < 1) %>% arrange(ct_FDR)
    print(paste("enriched RDI':", nrow(s1_rank), sep=' '))
  }
  if (direction == 'RUI'){
    print('RUI direction')
    s1_rank = fisher_output %>% subset(ct_FDR < 0.2 & OR > 1) %>% arrange(ct_FDR)
    print(paste("enriched RUI'", nrow(s1_rank), sep=' '))
  }
  
  return(s1_rank)
}

## step1_fisher(): function for differential activation analysis of (single) cci between on-treatment non-responder and pre-treatment samples (using two-sided fisher test)

  ## input:
    ## 1. activity_profile: binary vector of activated (1) or inactivated (0) cci; named vector by sample name
    ## 2. meta: data.frame with columns including group and sample; group corresponds to either sample belongs to on-treatment non-responder (1) or pre-treatment (0)

  ## output:
    ## data.frame with column p and OR; p corresponds to p-value; OR corresponds to odds-ratio; Odds-ratio > 1 suggest enrichment in on-treatment non-responder group

step1_fisher <- function(activity_profile, meta){ #function for fisher-exact test regarding enrichment of activation 
  
  on_nr = activity_profile[names(activity_profile) %in% subset(meta, group==1)$sample] #vector of all on-treatment non-responder (group = 1) samples (1: active; 0: inactive)
  pre = activity_profile[names(activity_profile) %in% subset(meta, group==0)$sample] #vector of all pre-treatment (group = 0) samples (1: active; 0: inactive)
  
  mat = data.frame(
    "present" = c(length(subset(on_nr, on_nr==1)), length(subset(pre, pre==1))), 
    "absent" = c(length(subset(on_nr, on_nr==0)), length(subset(pre, pre==0))),
    row.names= c("on_nr", "pre"),
    stringsAsFactors = F
  ) %>% as.matrix(.)
  # 4 cell matrix: 
  #       present absent
  # ON_NR O is 1  O is 0
  # PRE   P is 1  P is 0
  
  res = data.frame(p=fisher.test(mat)$p.value, OR = fisher.test(mat)$estimate) #fisher-exact test
  
  return(res) 
}

## ------- Functions for Step 2 -------
## step2(): function used to perform step 2 feature selection (greedy forward feature selection in pre-rx)
  ## input:
    ## s1_output: data.frame with columns including Interaction; Interaction contains the names of cci; data.frame should be pre arranged/sorted by p-value!!; Direct output of step 1 feature selection
    ## s2_ap: data.frame with col corresponding to cci name and rows [not named] corresponding to the same sample order in meta file (step 2 cohort activity profile)
    ## s2_meta: data.frame with columns including sample and response; sample corresponds to sample name ordered in parallel with the train_ap; response corresponds to ICB response: responder (1) and non-responder (0)
    ## trial: the number of times to iterate the greedy feature selection algorithm
    ## direction: direction to search for depleted ("RD") or activated ("RA") resistance relevant interactions

  ## output: 
    ## 2D matrix with rows as trials and columns as cci; value is feature score: 1 (selected by FFS and performed well on testing fold); 0 (not selected by FFS OR had "random" performance on testing fold); -1 (selected by FFS but performed poorly on testing fold)

step2 <- function(s1_output, s2_ap, s2_meta, trial=500, direction){
  
  s2_meta$response = factor(s2_meta$response, levels=c(0,1)) #recode response (level=1) and non-response (level=0) as factors
  auc_direction = ifelse(direction=='RDI', '<', '>') #set AUC direction based on whether we are looking for RDIs or RUIs

  fs_matrix = mclapply(1:trial, function(i){ #develop 2D feature score (0;1;-1) matrix iterating forward feature selection algorithm over n trials 
    
    print(i) # trial i
    
    #1. partition step 2 training cohort into 2 training and 1 testing fold 
    fold = KFold(c(1:nrow(s2_meta)), nfolds=3, seed=i)[[1]] 
    
    train_ap = s2_ap[-fold,]
    train_meta = s2_meta[-fold,]
    
    test_ap = s2_ap[fold,] 
    test_meta = s2_meta[fold,]
    
    if (length(unique(train_meta$response)) != 2 | length(unique(test_meta$response)) != 2){ #skip training with ffs algorithm if all training OR testing fold samples have the same response
      feature_score=rep(0, nrow(s1_output)) #feature score = 0 for all step 1 cci
      
      print('step2_train OR test have the same response information; cannot perform ffs; all cci recieved a feature score of 0 for this iteration') #warning message
      return(feature_score)
    }
    
    #2. run forward feature selection (FFS) algorithm on training fold
    train_index = step2_ffs(s1_output,train_ap,train_meta,auc_direction) #index of interactions selected from FFS on training fold 
    
    #3. evaluate performance on testing fold of interaction(s) selected from FFS using training fold
    test_df = test_ap %>% .[colnames(.) %in% s1_output[train_index,]$Interaction] %>% apply(.,1,sum) %>% data.frame(test_meta,sample_score=.) #calculate sample_score: sum of active queried interactions (from FFS) in each testing sample 
    res=suppressWarnings(suppressMessages(as.numeric(pROC::auc(as.numeric(test_df$response), as.numeric(test_df$sample_score), direction = auc_direction)))) #AUC (R vs. NR) using sum of active queried interactions (from FFS) in testing fold
    
    #4. assign feature score for selected interactions (from training fold) based on performance in testing fold
    feature_score=rep(0, nrow(s1_output)) #initialize all step 1 interactions (selected/not selected from FFS) w/ feature score of 0 
    
    if(res >= 0.6){ #if FFS selected interactions result in AUC >= 0.6, FFS selected interactions get a feature_score = 1
      feature_score[train_index] = 1 
    }
    if(res < 0.4){ #if ffs selected interactions result in AUC < 0.4, FFS selected interactions get a feature_score = -1
      feature_score[train_index] = -1
    }
    return(feature_score)
  })
  
  fs_matrix=do.call(rbind,fs_matrix)
  
  return(fs_matrix)
}

## step2_ffs(): forward feature selection algorithm to optimize power of cci in classifying responder and non-responder in training sample fold
  ## input:
    ## s1_output: data.frame with columns including Interaction (Lcell_Rcell_Ligand_Receptor); data.frame should be pre arranged (lowest ti greatest) by FDR per cell pair (ct_FDR)!!; Direct output of step 1 feature selection
    ## train_ap: data.frame with col corresponding to cci name and rows [not named] corresponding to the same sample order in meta file 
    ## train_meta: data.frame with columns including sample and response; sample corresponds to sample name ordered in parallel with the train_ap; response corresponds to ICB response: responder (1) and non-responder (0)
    ## auc_direction: direction to calculate auc; for resistant depleted "<" and for resistant activated ">"

  ## output: 
    ## numeric vector of indices regarding interactions selected by the greedy algorithm (selected_index)

step2_ffs <- function(s1_output, train_ap, train_meta, auc_direction){
  auc=0 #initalize start AUC
  current_index=0 #initialize index of the current feature being considered
  selected_index=c() #intializes indices of the selected features by the algorithm
  iteration = 0 #initialize counter for the number of iterations
  
  #forward greedy selection algorithm
  while((auc < 1) & iteration < nrow(s1_output)){ #iterate until AUC reaches 1 or maximum iterations reached (# interactions from step1)
    
    current_index = current_index + 1  #feature index
    query_cci = c(selected_index, current_index) %>% s1_output[.,] #select features based on prev. selected indices by the algorithm and current index
    train_df = train_ap %>% .[colnames(.) %in% query_cci$Interaction] %>% apply(.,1,sum) %>% data.frame(train_meta,sample_score=.) #calculate sample_score: sum of active queried interactions in each training sample 
    res=suppressWarnings(suppressMessages(as.numeric(pROC::auc(as.numeric(train_df$response), as.numeric(train_df$sample_score), direction = auc_direction)))) #AUC (R vs. NR) using sum of active queried interactions in training fold
    
    #evaluate performance of queried interactions 
    if (res > auc){ #if AUC from queried features better than selected indices 
      selected_index = c(selected_index, current_index) #append current feature (index) to selected indices by algorithm
      auc = res #update best AUC 
    }
    
    iteration = iteration + 1
    
  }
  return(selected_index)  #return the indices of the selected features
}

## ------- Additional Functions -------

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## ------- START CODE HERE -------

## parallelize code
ncpus <- parallelly::availableCores()
print(ncpus)
options(mc.cores = ncpus) # set a global option for parallel packages

## training input
model_input = loadRData("/IRIS data/LIRICS/organized data/ICB/IRIS_LIRICS_final_ICB_cohort_input_response.Rdata") #CHANGE
test_input = model_input[[1]] #test input(s)
s1_input = model_input[[1]] #step 1 input(s)
s2_input = model_input[[1]][c(2,3,5,8)] #step 2 input(s)

## Perform LOOCohort Validation with IRIS
direction='RDI' #RDI for resistance downregulated; RUI for Resistance Upregulated
ensembl=lapply(1:length(test_input), function(e) LOOC_validation(names(test_input[e]), s1_input, s2_input, direction=direction))
names(ensembl) = names(test_input)

## save file
save(file="XXX", ensembl)
