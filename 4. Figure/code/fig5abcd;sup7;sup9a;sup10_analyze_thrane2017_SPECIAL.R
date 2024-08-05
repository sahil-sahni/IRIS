library(tidyverse)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(pROC)
library(rstatix)
library(scales)
library(abind)
library(pROC)
library(viridis)


## ------- Function -------
## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

#read input Thrane et al. 2018 file
input = thrane_output

# load benchmark signature
signature = sig_thrane

# load RDI interactions
RDI = feature_RDI
RDI = RDI %>% unlist(.) %>% unique(.)

## ------- Start -------

# organize data
RDS_df = lapply(1:length(input), function(x){
  
  slidescore = lapply(1:length(input[[x]]), function(y){ #length(input[[x]])
    
    # load slide information
    name = input[[x]][[y]][[1]]
    print(name)
    
    name2 = substr(name, 1, nchar(name) - 7)
    
    if(name %in% names(signature)){
    icb = signature[[name]] %>% mutate(Row.names=gsub("X", "", Row.names))
    }
    
    if(!(name %in% names(signature))){
      icb = signature[[name2]] %>% mutate(Row.names=gsub("X", "", Row.names))
    }
    
    cf = input[[x]][[y]][[2]]
    
    if (!('TCD8' %in% colnames(cf))){
      cf = cf %>% mutate(TCD8 = 0)
    }
    
    cf = cf %>% select(TCD8) %>% mutate(TCD8_AUC=ifelse(TCD8 == 0, 0, 1)) #cell fraction 
    
    # run SOCIAL/SPECIAL step 4
    spotscore = SOCIAL.infer_activity(input[[x]][[y]]) %>% t(.)

    # subset for RDI relevant interactions
    spotscore_RDI = spotscore %>% subset(rownames(.) %in% RDI)

    # calculate RDS score and create RDS df with all individual interactions
    RDI_score <- spotscore_RDI %>% colSums(.,na.rm=T) %>% as.data.frame(.) %>% set_colnames('RDI')
    output = merge(cf,RDI_score, by.x=0, by.y=0) %>% set_colnames(c('coord','TCD8 fraction','TCD8 infiltration','active_RDI'))
    
    # for plotting RDS score
    output = output %>% select(1:4) %>%
      mutate(
        sample = x,
        name = input[[x]][[y]][[1]],
        slide = ifelse(grepl('null',name)==F,'foreground','background') %>% factor(., levels=c('foreground','background')),
        x = as.numeric(str_split(coord, "x") %>% sapply(function(x) x[1])),
        y = as.numeric(str_split(coord, "x") %>% sapply(function(x) x[2])),
        `TCD8 infiltration` = factor(as.character(`TCD8 infiltration`), levels=c(1,0)),
        RDS_scale = c(scale(active_RDI)),
        `TCD8 fraction_scale` = c(scale(`TCD8 fraction`))
      ) 
    
    output = merge(output, icb, by.x=1, by.y=1)

    
    

    print(cor(output$active_RDI, output$`TCD8 fraction`, method = c('spearman')))
    return(output)
  })
  
  res = do.call(rbind, slidescore)
  return(res)
})
  
res = do.call(rbind,RDS_df) #%>% subset(slide != 'background') # remove background slides for all
res_tcd8 = res %>% subset(sample %in% c(1,2,4,5)) #slides with TCD8 infiltratoin
res_tcd8_v2 = res_tcd8 %>% mutate(patient = substr(name, 1, nchar(name) - 5)) #patient level for slides with TCD8 infiltration
res_woutb = res %>% mutate(group=ifelse(sample %in% c(1,2,4,5), 'w/CD8+', 'w/out CD8+')) #all slides

## ------- Figure 5C per patient level Thrane et al. 2018 -------

stat.test = res_tcd8_v2 %>% mutate(ti=`TCD8 infiltration`) %>% group_by(patient) %>% wilcox_test(active_RDI ~ ti, ref.group = '1', alternative = 'greater')
stat.test <- stat.test %>%
  add_xy_position(x = "patient", dodge = 0.8) %>% add_significance(p.col='p')
fig5c = ggboxplot(res_tcd8_v2 %>% mutate(ti=`TCD8 infiltration`, ti_binary = ifelse(ti == 1,'yes','no')) %>% mutate(ti_binary=factor(ti_binary,levels=c('yes','no'))), x='patient', y='active_RDI', fill='ti_binary', size=0.3, outlier.size=1) + 
  theme(legend.position = 'right', axis.title.x = element_blank(), axis.text.x = element_text(angle=-30,hjust=0))
fig5c = fig5c + stat_pvalue_manual(stat.test,  label = "p", tip.length = .02, y.position=17, bracket.size = 0.3, size=10/.pt)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels=label_number(accuracy=1))+
  scale_fill_manual(values=c('#D0759F','#88BECD'))+
  ylab('# of active RDIs')+
  labs(fill = "CD8+ T cell infiltration")
plot(fig5c)

## ------- Supp 9A per slide level Thrane et al. 2018 -------

stat.test = res_tcd8_v2 %>% mutate(ti=`TCD8 infiltration`) %>% group_by(name) %>% wilcox_test(active_RDI ~ ti, ref.group = '1', alternative = 'greater')
stat.test <- stat.test %>%
  add_xy_position(x = "name", dodge = 0.8) %>% add_significance(p.col='p')
sup9a = ggboxplot(res_tcd8_v2 %>% mutate(ti=`TCD8 infiltration`, ti_binary = ifelse(ti == 1,'yes','no')) %>% mutate(ti_binary=factor(ti_binary,levels=c('yes','no'))), x='name', y='active_RDI', fill='ti_binary', size=0.3, outlier.size=1) + 
  theme(legend.position = 'right', axis.title.x = element_blank(), axis.text.x = element_text(angle=-30,hjust=0))
sup9a = sup9a + stat_pvalue_manual(stat.test,  label = "p", tip.length = .02, y.position=17, bracket.size = 0.3, size=10/.pt)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels=label_number(accuracy=1))+
  scale_fill_manual(values=c('#D0759F','#88BECD'))+
  ylab('# of active RDIs')+
  labs(fill = "CD8+ T cell infiltration")
plot(sup9a)

# ------- Figure 5D comparing slides w/ and w/out CD8+ T cell infiltration -------
stat.test = res_woutb  %>% wilcox_test(active_RDI ~ group, ref.group = 'w/CD8+', alternative = 'greater')
stat.test <- stat.test %>%
  add_xy_position(x = "group", dodge = 0.8) %>% add_significance(p.col='p')
fig5d = ggviolin(res_woutb, x='group', y='active_RDI', fill='group', size=0.3) + 
  theme(legend.position = 'none', axis.title.x = element_blank(), axis.text.x = element_text(angle=-30,hjust=0))
fig5d = fig5d + stat_pvalue_manual(stat.test,  label = "p", tip.length = .02, y.position=18, bracket.size = 0.3, size=10/.pt)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels=label_number(accuracy=1))+
  scale_fill_manual(values=c('#D0759F','#88BECD'))+
  ylab('# of active RDIs')+
  labs(fill = "CD8+ T cell infiltration")
plot(fig5d)

## ------- Figure 5A Spatial Map of RDS interaction -------
fig5a = ggplot(res_tcd8, aes(x,y))+
  geom_point(aes(color=RDS_scale), size=1)+
  scale_color_viridis(option='magma', direction=-1)+
  facet_wrap(~name, nrow=2, ncol=2, strip.position = c('top'))+
  theme_pubr()+
  labs(color= '# of active RDIs (scaled)')
plot(fig5a)

## ------- Figure 5B Spatial Map of RDS interaction -------

fig5b = ggplot(res_tcd8, aes(x,y))+
  geom_point(aes(color=`TCD8 fraction_scale`), size=1)+
  scale_color_viridis(option='magma', direction=-1)+
  facet_wrap(~name, nrow=2, ncol=2)+
  theme_pubr()+
  labs(color= 'fraction of CD8+ T cells (scaled)')
plot(fig5b)

## ------- Sup 10A Comparison of RDI activity with ICB signature -------

# calculate signature score
res_sig_v3 = res %>% mutate(RDI_present=ifelse(active_RDI > 1, 'yes','no')) %>% pivot_longer(c(16:21), names_to='signature', values_to='score') %>%
  mutate(score = ifelse(signature == 'MPS', score * -1, score)) %>% group_by(signature) %>% mutate(score_scale = c(scale(score))) %>% ungroup() %>% group_by(sample,coord) %>% mutate(average_score=mean(score_scale) ) %>%
  mutate(RDI_present=factor(RDI_present, levels=c('yes','no')), signature=factor(signature))%>% mutate(patient = substr(name, 1, nchar(name) - 5)) 


stat.test = res_sig_v3 %>% group_by(patient,signature) %>% wilcox_test(score ~ RDI_present, ref.group = 'yes', alternative = 'greater')
stat.test <- stat.test %>%
  add_xy_position(x = "patient", dodge = 0.8) %>% add_significance(p.col='p')
sup10a = ggboxplot(res_sig_v3 , x='patient', y='score', fill='RDI_present', facet.by = 'signature', size=0.3, outlier.size=1) + 
  theme(legend.position = 'bottom', axis.title.x = element_blank(), axis.text.x = element_text(angle=-30,hjust=0))+
  facet_wrap(~signature, scales='free')
sup10a = sup10a + stat_pvalue_manual(stat.test,  label = "p", tip.length = .02, bracket.size = 0.3, size=7/.pt, y.position = c(1,3,6.5,7.4,2500,2.2))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels=label_number(accuracy=.1))+
  scale_fill_manual(values=c('#D0759F','#88BECD'))+
  ylab('ICB signature score')+
  labs(fill = "RDI present")
plot(sup10a)

## ------- Sup 7A Spatial Map of RDS interaction (all slides) -------
sup7a = ggplot(res, aes(x,y))+
  geom_point(aes(color=RDS_scale), size=1)+
  scale_color_viridis(option='magma', direction=-1)+
  facet_wrap(~name, nrow=2, ncol=4, strip.position = c('top'))+
  theme_pubr()+
  labs(color= '# of active RDIs (scaled)')
plot(sup7a)

## ------ Calculate AUC to test overlap of pattern ------

df = res_tcd8 %>% mutate(ti = factor(`TCD8 infiltration`, levels=c(0,1)))

# Create an empty vector to store AUC values
auc_values <- numeric()
# Loop through each unique sample
for (s in unique(df$sample)) {
  # Subset data for the current sample
  data_subset <- df[df$sample == s, ]
  
  # Calculate AUC curve
  auc_value <- pROC::auc(data_subset$ti, data_subset$active_RDI, direction='<')
  # Print or store the AUC value
  print(paste("Sample", s, "AUC:", auc_value))
  
  # Store the AUC value in the vector
  auc_values <- c(auc_values, auc_value)
}

# Print AUC values for all samples
print("AUC values for all samples:")
print(auc_values)
median(auc_values)
sd(auc_values)


