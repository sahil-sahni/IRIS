Code to reproduce the results of existing biomarkers of ICB response. The TPM calculation for Auslander et al. and PUCH were used for both ICB tool benchmarking and CODEFACS deconvolution.

**Auslander et al.**:
To reproduce the bulk transcriptomics results for Auslander et al., counts data was first converted to TPM following this code:

```r
library(tidyverse)

#Auslander et al. 2018
a2018 = read.csv('GSE115821_MGH_counts.csv', check.names = F) #from GEO GSE115821
a2018_count = a2018[,7:43]
a2018_length = a2018[,6]
a2018_gene = a2018 %>% dplyr::select(Gene = Geneid)
x = a2018_count/a2018_length 
a2018_tpm <- t( t(x) * 1e6 / colSums(x) )
a2018_final = cbind(a2018_gene, a2018_tpm)
```

**PUCH**:
To reproduce the bulk transcriptomics results for PUCH, FPKM data was first converted to TPM following this code:

```r
#PUCH
library(parallel)
library(tidyverse)

puch <- read.csv("melanoma_puch_exp.csv") # from Cui et al. 
puch_gene = puch %>% dplyr::select(Gene=X)
puch_fpkm = (2^(puch[,-1]))-1
puch_tpm = mclapply(1:nrow(puch_fpkm), function(x){
  print(x)
  fpkm_i = puch_fpkm[x,]
  fpkm_j = colSums(puch_fpkm)
  tpm = (fpkm_i/fpkm_j)*10^6
})
puch_tpm = do.call(rbind, puch_tpm)
colSums(puch_tpm)
puch_tpm = cbind(puch_gene, puch_tpm)
```

**Thrane et al.**:
To reproduce the bulk transcriptomics results for Thrane et al., counts data was first converted to TPM using in-house code. In brief, counts data was converted to TPM using gene lengths from GRCh37.87. Duplicate genes and non-protein coding genes were removed, and data was renormalized to TPM with unique pcg.
