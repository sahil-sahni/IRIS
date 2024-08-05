calculate_score <- function(test, name, features){ #
  input_name = strsplit(name, '_') %>% unlist(.) %>% .[[1]] 
  features_name = strsplit(names(features),'_') %>% sapply(.,'[[',1) %>% unlist(.)
  
  train = features[which(input_name != features_name)] #ensure ensembl model is not training and testing with the same cohort
  
  ensembl_df=lapply(1:length(train), function(r){ #calculate score for each ensembl model
    model_score = test[,-1:-2] %>% .[colnames(.) %in% train[[r]]] %>% apply(1, sum) 
    model_score = model_score/length(train[[r]])
    return(model_score)
  })
  ensembl_df=do.call(cbind, ensembl_df)
  
  score_df = ensembl_df %>% apply(1,function(x){sum(x,na.rm=T)}) %>% data.frame(response=test[,2], score=./length(train)) ###NEW /length(train)
  
  return(score_df)
}
calculate_score_surv <- function(test, name, features){ #
  input_name = strsplit(name, '_') %>% unlist(.) %>% .[[1]] 
  features_name = strsplit(names(features),'_') %>% sapply(.,'[[',1) %>% unlist(.)
  
  train = features[which(input_name != features_name)] #ensure ensembl model is not training and testing with the same cohort
  
  ensembl_df=lapply(1:length(train), function(r){ #calculate score for each ensembl model
    model_score = test[,-1:-4] %>% .[colnames(.) %in% train[[r]]] %>% apply(1, sum) 
    model_score = model_score/length(train[[r]])
    return(model_score)
  })
  ensembl_df=do.call(cbind, ensembl_df)
  
  score_df = ensembl_df %>% apply(1,function(x){sum(x,na.rm=T)}) %>% data.frame(test[,1:4], score=./length(train))
  
  return(score_df)
}
calculate_score_scaled <- function(test, name, features){
  input_name = strsplit(name, '_') %>% unlist(.) %>% .[[1]] 
  features_name = strsplit(names(features),'_') %>% sapply(.,'[[',1) %>% unlist(.)
  
  train = features[which(input_name != features_name)] #ensure ensembl model is not training and testing with the same cohort
  
  ensembl_df=lapply(1:length(train), function(r){ #calculate score for each ensembl model
    model_score = test[,-1:-2] %>% .[colnames(.) %in% train[[r]]] %>% apply(1, sum) 
    model_score = model_score/length(train[[r]])
    return(model_score)
  })
  ensembl_df=do.call(cbind, ensembl_df)
  
  score_df = ensembl_df %>% apply(1,function(x){sum(x,na.rm=T)}) %>% data.frame(response=test[,2], score=./length(train), scale=scale(./length(train)))
  
  return(score_df)
}
calculate_score_scaled_tcga <- function(test, name, features){
  input_name = strsplit(name, '_') %>% unlist(.) %>% .[[1]] 
  features_name = strsplit(names(features),'_') %>% sapply(.,'[[',1) %>% unlist(.)
  
  train = features[which(input_name != features_name)] #ensure ensembl model is not training and testing with the same cohort
  
  ensembl_df=lapply(1:length(train), function(r){ #calculate score for each ensembl model
    model_score = test[,-1:-16] %>% .[colnames(.) %in% train[[r]]] %>% apply(1, sum) 
    model_score = model_score/length(train[[r]])
    return(model_score)
  })
  ensembl_df=do.call(cbind, ensembl_df)
  
  score_df = ensembl_df %>% apply(1,function(x){sum(x,na.rm=T)}) %>% data.frame(test[,1:16], score=., scale=scale(./length(train)))
  
  return(score_df)
}
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
