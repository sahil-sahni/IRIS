library(tidyverse)
library(pROC)
library(magrittr)
library(ggpubr)
library(rlist)
library(arules)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# function to calculate AUC based on signature score
calc_auc <- function(score, cohort_name, meta, direction='<'){
  meta = meta[[cohort_name]] %>% dplyr::select(sample,response) %>% merge(., as.data.frame(score), by.x=1, by.y=0)

  #calculate auc
  auc = suppressWarnings(suppressMessages(as.numeric(pROC::auc(as.numeric(meta$response), as.numeric(meta$score), direction = direction))))
  
  return(auc)
}

## organize bulk dataset

#spatial data
m1r1 = read_tsv('ST_mel1_rep1_TPM.tsv') %>% as.data.frame()
m1r1 = m1r1 %>% set_rownames(., .$Gene) %>% .[,-1] %>% as.matrix(.)

m1r2 = read_tsv('ST_mel1_rep2_TPM.tsv') %>% as.data.frame()
m1r2 = m1r2 %>% set_rownames(., .$Gene) %>% .[,-1] %>% as.matrix(.)

m2r1 = read_tsv('ST_mel2_rep1_TPM.tsv') %>% as.data.frame()
m2r1 = m2r1 %>% set_rownames(., .$Gene) %>% .[,-1] %>% as.matrix(.)

m2r2 = read_tsv('ST_mel2_rep2_TPM.tsv') %>% as.data.frame()
m2r2 = m2r2 %>% set_rownames(., .$Gene) %>% .[,-1] %>% as.matrix(.)

m3r1 = read_tsv('ST_mel3_rep1_TPM.tsv') %>% as.data.frame()
m3r1 = m3r1 %>% set_rownames(., .$Gene) %>% .[,-1] %>% as.matrix(.)

m3r2 = read_tsv('ST_mel3_rep2_TPM.tsv') %>% as.data.frame()
m3r2 = m3r2 %>% set_rownames(., .$Gene) %>% .[,-1] %>% as.matrix(.)

m4r1 = read_tsv('ST_mel4_rep1_TPM.tsv') %>% as.data.frame()
m4r1 = m4r1 %>% set_rownames(., .$Gene) %>% .[,-1] %>% as.matrix(.)

m4r2 = read_tsv('ST_mel4_rep2_TPM.tsv') %>% as.data.frame()
m4r2 = m4r2 %>% set_rownames(., .$Gene) %>% .[,-1] %>% as.matrix(.)

bulk_datasets = list(mel1_rep1=m1r1, mel1_rep2=m1r2, mel2_rep1=m2r1, mel2_rep2=m2r2, mel3_rep1=m3r1, mel3_rep2=m3r2, mel4_rep1=m4r1, mel4_rep2=m4r2)

##################################################################################################
#Calculate gene signature scores (PD1, PDL1, CTLA4, Texhaustion)

gene_score <- function(bulk){
  bulk <- as.data.frame(bulk)
  gene_sig <- list(PD1=c('PDCD1'), PDL1=c('CD274'), CTLA4=c('CTLA4'), Tex = c('LAG3','TIGIT','PDCD1','PD1','HAVCR2','TIM3','CTLA4'))
  score_list = list()

  for (i in 1:length(gene_sig)){
    
    gene_sig_fil = gene_sig[[i]] %>% subset(. %in% row.names(bulk))
    s1 <- colMeans(bulk[gene_sig_fil,],na.rm = F)
    
    #zscore computation
    scores = sapply(s1, function(x) (x - mean(s1,na.rm = T))/sd(s1, na.rm = T))
    
    score_list = list.append(score_list, s1)
    names(score_list)[i] = names(gene_sig[i])
  }
  
  return(score_list)
}
# calculate signature score for each bulk data set
scores_gene = sapply(bulk_datasets, function(x) gene_score(x))


##################################################################################################
#MPS score calculator Merlino et al 2020
MPS_score <- function(bulk){
  bulk <- as.data.frame(bulk)
  MPS_plus <- c("AKR1C3","BMP1","CRTAC1","ECEL1","ERC2","FAM110C","FUT9","GABRA2","GAP43","GREM1","HECW1","KLHL1","KRT12","LHFPL4","NEFL","NEFM","NETO1","NKX2-2","NSG2","OCIAD2","OTOP1","PDE3B","PTPRN2","PTPRT","SIGLEC15","SLC13A5","SLC9A2","SLITRK6","SNAP91","STON2","TAC1","VAT1L","WNT5A")
  MPS_minus <- c("ALX1","BRD7","DTD1","GRSF1","HCN1","LTA4H","OXCT1","PATJ","PLXNC1","SSBP4","TELO2","TMEM177")
  s1 <- colSums(bulk[MPS_plus,],na.rm = T)
  s2 <- colSums(bulk[MPS_minus,], na.rm = T)
  
  #zscore computation
  scores = sapply(s1-s2, function(x) (x - mean(s1-s2,na.rm = T))/sd(s1-s2, na.rm = T))
  return(scores)
}

# calculate MPS score
scores_MPS <- sapply(bulk_datasets, function(x) MPS_score(x))


##################################################################################################
#IFNG signature score calculator Fehrenbacher et al 2016, Ayers 2017
IFNG_sig_score <- function(bulk){
  bulk <- as.data.frame(bulk)
  IFNG_sig <- c('IFNG', 'STAT1', 'HLA-DRA', 'CXCL9', 'CXCL10', 'IDO1')
  house_keeping <- c("STK11IP","ZBTB34","TBC1D10B","OAZ1","POLR2A","G6PD","ABCF1","C14orf102","UBB","TBP","SDHA")
  s1 <- log10(bulk[IFNG_sig,]+1)
  s2 <- log10(bulk[house_keeping,]+1)
  s1 <- s1 - rowMeans(s2, na.rm = T)
  
  #score computation
  scores = colMeans(s1, na.rm = T)
  return(scores)
}

# calculate IFN score
scores_IFNG <- sapply(bulk_datasets, function(x) IFNG_sig_score(x))

##################################################################################################
#Davoli cytolytic gene signature score

Davoli_score <- function(bulk){
  
  # Literature genes
  Davoli_IS.read <- c("CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK")
  match_Davoli_IS.genes <- match(Davoli_IS.read, rownames(bulk))
  
  # Log2 transformation:
  log2.RNA.tpm <- log2(bulk + 1)
  
  # Subset log2.RNA.tpm
  sub_log2.RNA.tpm <- log2.RNA.tpm[match_Davoli_IS.genes, ]
  
  # Calculate rank position for each gene across samples
  ranks_sub_log2.RNA.tpm <- apply(sub_log2.RNA.tpm, 1, rank)
  
  # Get normalized rank by divided
  ranks_sub_log2.RNA.tpm.norm <- (ranks_sub_log2.RNA.tpm - 1)/(nrow(ranks_sub_log2.RNA.tpm) - 1)
  
  # Calculation: average of the expression value of all the genes within-sample
  score <- apply(ranks_sub_log2.RNA.tpm.norm, 1, mean)
  return(score)
}

# calculate Davoili score
scores_Davoli <- sapply(bulk_datasets, function(x) Davoli_score(x))

##################################################################################################
#T cell inflamed GEP
Tcell_INF_score <- function(bulk){
  bulk <- as.data.frame(bulk)
  Tcell_inf_sig <- c("CXCR6","TIGIT","CD27","CD274","PDCD1LG2","LAG3","NKG7","PSMB10","CMKLR1","CD8A","IDO1","CCL5","CXCL9","HLA-DQA1","CD276", "HLA-DRB1","HLA-E","STAT1")
  house_keeping <- c("STK11IP","ZBTB34","TBC1D10B","OAZ1","POLR2A","G6PD","ABCF1","C14orf102","UBB","TBP","SDHA")
  s1 <- log10(bulk[Tcell_inf_sig,]+1)
  s2 <- log10(bulk[house_keeping,]+1)
  s1 <- s1 - rowMeans(s2, na.rm = T)
  
  weights = c(0.004313, 0.084767, 0.072293, 0.042853, 0.003734, 0.123895, 0.075524, 0.032999, 0.151253, 0.031021, 0.060679, 0.008346, 0.074135, 0.020091, 0.058806, 0.07175, 0.250229)
  s1 <- s1*weights
  scores <- colSums(s1, na.rm = T)
  return(scores)
  
}

# calculate inflamed GEP score
scores_Tin <- sapply(bulk_datasets, function(x) Tcell_INF_score(x))


##################################################################################################
#OE resistance program from Livnat et al. 2018
discretize<-function(v,n.cat){
  q1<-quantile(v,seq(from = (1/n.cat),to = 1,by = (1/n.cat)))
  u<-matrix(nrow = length(v))
  for(i in 2:n.cat){
    u[(v>=q1[i-1])&(v<q1[i])]<-i
  }
  return(u)
}
add.up.down.suffix<-function(v){
  v<-lapply(v, function(x) x<-c(x,paste0(x,".up"),paste0(x,".down")))
  v<-unlist(v,use.names = F)
  return(v)
}
get.OE.bulk <- function(r,gene.sign = NULL,num.rounds = 1000,full.flag = F){
  set.seed(1234)
  r$genes.mean<-rowMeans(r$tpm)
  r$zscores<-sweep(r$tpm,1,r$genes.mean,FUN = '-')
  r$genes.dist<-r$genes.mean
  r$genes.dist.q<-discretize(r$genes.dist,n.cat = 50)
  r$sig.scores<-matrix(data = 0,nrow = ncol(r$tpm),ncol = length(gene.sign))
  sig.names<-names(gene.sign)
  colnames(r$sig.scores)<-sig.names
  r$sig.scores.raw<-r$sig.scores
  rand.flag<-is.null(r$rand.scores)|!all(is.element(names(gene.sign),colnames(r$rand.scores)))
  if(rand.flag){
    print("Computing also random scores.")
    r$rand.scores<-r$sig.scores
  }
  for (i in sig.names){
    b.sign<-is.element(r$genes,gene.sign[[i]])
    if(sum(b.sign)<2){next()}
    if(rand.flag){
      rand.scores<-get.semi.random.OE(r,r$genes.dist.q,b.sign,num.rounds = num.rounds)
    }else{
      rand.scores<-r$rand.scores[,i]
    }
    raw.scores<-colMeans(r$zscores[b.sign,])
    final.scores<-raw.scores-rand.scores
    r$sig.scores[,i]<-final.scores
    r$sig.scores.raw[,i]<-raw.scores
    r$rand.scores[,i]<-rand.scores
  }
  if(full.flag){return(r)}
  sig.scores<-r$sig.scores
  #sig.scores = sig.scores[,"fnc.up"] - sig.scores[,"fnc.down"]
  #names(sig.scores) = r$samples
  return(sig.scores)
}
get.semi.random.OE <- function(r,genes.dist.q,b.sign,num.rounds = 1000,full.flag = F){
  # Previous name: get.random.sig.scores
  sign.q<-as.matrix(table(genes.dist.q[b.sign]))
  q<-rownames(sign.q)
  idx.all<-c()
  B<-matrix(data = F,nrow = length(genes.dist.q),ncol = num.rounds)
  Q<-matrix(data = 0,nrow = length(genes.dist.q),ncol = num.rounds)
  for (i in 1:nrow(sign.q)){
    num.genes<-sign.q[i]
    if(num.genes>0){
      idx<-which(is.element(genes.dist.q,q[i]))
      for (j in 1:num.rounds){
        idxj<-sample(idx,num.genes) 
        Q[i,j]<-sum(B[idxj,j]==T)
        B[idxj,j]<-T
      }  
    }
  }
  rand.scores<-apply(B,2,function(x) colMeans(r$zscores[x,]))
  if(full.flag){return(rand.scores)}
  rand.scores<-rowMeans(rand.scores)
  return(rand.scores)
}
cmb.res.scores<-function(r,res.sig = NULL,bulk.flag = T){
  res<-r$res
  if(!is.null(res.sig)){
    res<-get.OE(r = r,sig = res.sig,bulk.flag = bulk.flag)
  }
  
  if(!all(is.element(c("trt","exc","exc.seed"),colnames(res)))){
    print("No merge.")
    return(res)
  }
  res<-cbind.data.frame(excF = rowMeans(res[,c("exc","exc.seed")]),
                        excF.up = rowMeans(res[,c("exc.up","exc.seed.up")]),
                        excF.down = rowMeans(res[,c("exc.down","exc.seed.down")]),
                        res = rowMeans(res[,c("trt","exc","exc.seed")]),
                        res.up = rowMeans(res[,c("trt.up","exc.up","exc.seed.up")]),
                        res.down = rowMeans(res[,c("trt.down","exc.down","exc.seed.down")]),
                        res)
  
  if(is.element("fnc",colnames(res))){
    res<-cbind.data.frame(resF = res[,"res"]+res[,"fnc"],
                          resF.up = res[,"res.up"]+res[,"fnc.up"],
                          resF.down = res[,"res.down"]+res[,"fnc.down"],
                          res)
  }
  
  remove.sig<-add.up.down.suffix(c("exc","exc.seed","fnc"))
  res<-res[,setdiff(colnames(res),remove.sig)]
  colnames(res)<-gsub("excF","exc",colnames(res))
  return(res)
}
compute.samples.res.scores<-function(r,res.sig=NULL,cell.sig=NULL,residu.flag = F,
                                     cc.sig = NULL,num.rounds = 1000){
  # modified input 
  res.sig=loadRData("~/Documents/anti_pdl1_project/5. benchmark/signatures/livnat jerby signatures/resistance.program.RData") # download data from https://github.com/livnatje/ImmuneResistance/blob/master/Results/Signatures/resistance.program.RData
  cell.sig=loadRData("~/Documents/anti_pdl1_project/5. benchmark/signatures/livnat jerby signatures/cell.type.sig.RData") # download data from https://github.com/livnatje/ImmuneResistance/blob/master/Results/Signatures/cell.type.sig.RData
  r=list(tpm=r, genes=row.names(r), samples=colnames(r))
  
  r$res.ori<-NULL;r$res<-NULL;r$tme<-NULL;r$X<-NULL
  names(res.sig)<-gsub(" ",".",names(res.sig))
  two.sided<-unique(gsub(".up","",gsub(".down","",names(res.sig))))
  b<-is.element(paste0(two.sided,".up"),names(res.sig))&
    is.element(paste0(two.sided,".down"),names(res.sig))
  two.sided<-two.sided[b]
  r$res<-get.OE.bulk(r,res.sig,num.rounds = num.rounds)
  res2<-r$res[,paste0(two.sided,".up")]-r$res[,paste0(two.sided,".down")]
  if(!is.matrix(res2)){res2<-as.matrix(res2)}
  colnames(res2)<-two.sided
  r$res<-cbind(res2,r$res)
  rownames(r$res)<-colnames(r$tpm)
  r$tme<-get.OE.bulk(r,cell.sig,num.rounds = num.rounds)
  r$res<-cmb.res.scores(r)
  if(residu.flag){
    print("Computing residuales")
    r$res<-t(get.residuals(t(r$res),colSums(r$tpm>0)))
    r$tme<-t(get.residuals(t(r$tme),colSums(r$tpm>0)))
  }
  if(!is.null(cc.sig)){
    print("Filtering cell-cycle effects")
    if(is.null(r$cc.scores)){
      r$cc.scores<-get.OE.bulk(r,cc.sig,num.rounds = num.rounds)
    }
    r$res.ori<-r$res
    r$res<-t(get.residuals(t(r$res),r$cc.scores))
  }
  if(is.element("resF",colnames(r$res))&is.element("T.CD8",colnames(r$tme))){
    r$res<-cbind.data.frame(r$res,resF.minus.TIL = r$res[,"resF"]-r$tme[,"T.CD8"])
    if(!is.null(r$res.ori)){
      r$res.ori<-cbind.data.frame(r$res.ori,resF.minus.TIL = r$res.ori[,"resF"]-r$tme[,"T.CD8"])
    }
  }
  
  ## modified returned output!
  score = r$res$resF
  names(score) = r$samples
  return(score)
}

# calculate functional resistance score
scores_resF <- sapply(bulk_datasets, function(x) compute.samples.res.scores(x))

##################################################################################################
#Testing predictive performance of IMPRES, Auslander et al 2018
#IMPRES score calculator
IMPRES_score <- function(bulk){
  bulk <- as.data.frame(bulk)
  impres_pairs <- matrix(c("PDCD1","TNFSF4",
                           "CD27","PDCD1",
                           "CTLA4","TNFSF4",
                           "CD40","CD28",
                           "CD86","TNFSF4",
                           "CD28","CD86",
                           "CD80","TNFSF9",
                           "CD274","VSIR",
                           "CD86","HAVCR2",
                           "CD40","PDCD1",
                           "CD86","CD200",
                           "CD40","CD80",
                           "CD28","CD276",
                           "CD40","CD274",
                           "TNFRSF14","CD86"),15,2)
  scores = rowSums(apply(impres_pairs,1,function(x) return(bulk[x[1],] > bulk[x[2],])),na.rm = T)
  
  ## modified returned output!
  names(scores)=colnames(bulk)
  return(scores)
}

# calculate IMPRES score
scores_IMPRES <- sapply(bulk_datasets, function(x) IMPRES_score(x))

# merge score input
merge_score = function(x,y,name,i,j){
  y=as.data.frame(y) %>% set_colnames(name)
  z= merge(x, y, by.x=i, by.y=j)
}

#for score analysis
score_list = scores_gene %>% apply(2, function(x){as.data.frame(x)}) %>%
  mapply(merge_score, ., scores_Davoli, MoreArgs=list(name='Cytotoxic',i=0,j=0)) %>% apply(2, function(x){as.data.frame(x)}) %>% 
  mapply(merge_score, ., scores_IFNG, MoreArgs=list(name='IFNG',i=1,j=0)) %>% apply(2, function(x){as.data.frame(x)}) %>% 
  mapply(merge_score, ., scores_MPS, MoreArgs=list(name='MPS',i=1,j=0)) %>% apply(2, function(x){as.data.frame(x)}) %>% 
  mapply(merge_score, ., scores_Tin, MoreArgs=list(name='Tin',i=1,j=0)) %>% apply(2, function(x){as.data.frame(x)}) %>%
  mapply(merge_score, ., scores_resF, MoreArgs=list(name='resF',i=1,j=0)) %>% apply(2, function(x){as.data.frame(x)})%>%
  mapply(merge_score, ., scores_IMPRES, MoreArgs=list(name='IMPRES',i=1,j=0)) %>% apply(2, function(x){as.data.frame(x)})


saveRDS(score_list, file='sig_input_spatial.RDS')


