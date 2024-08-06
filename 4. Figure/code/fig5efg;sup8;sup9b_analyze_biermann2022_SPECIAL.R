library(tidyverse)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(pROC)
library(rstatix)
library(scales)
library(abind)
library(viridis)

## ------- Function -------
## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

#read input biermann et al. 2022 file
input = biermann_output

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
    
    print(cor(output$active_RDI, output$`TCD8 fraction`, method = c('spearman')))
    
    return(output)
  })
  
  res = do.call(rbind, slidescore)
  return(res)
})
  
res = do.call(rbind,RDS_df)

#res_tcd8 = subset(res, slide != 'background') 
res = res %>% mutate(patient = substr(name, 1, nchar(name) - 2)) %>% subset(!(name %in% c('MBM13.1','MBM07.1','ECM08.1'))) %>% mutate(group=ifelse(patient %in% c('XX','XX'),'PD','naïve')) #XX information was retrieved from Ben Izar's Group

## ------- Figure 5E Boxplot score between spatial spots with and without cd8 infiltration patient level -------
stat.test = res %>% mutate(ti=`TCD8 infiltration`) %>% group_by(patient) %>% rstatix::wilcox_test(active_RDI ~ ti, ref.group = '1', alternative = 'greater')
stat.test <- stat.test %>%
  add_xy_position(x = "patient", dodge = 0.8) %>% add_significance(p.col='p')
fig5e = ggboxplot(res%>% mutate(ti=`TCD8 infiltration`, ti_binary = ifelse(ti == 1,'yes','no')) %>% mutate(ti_binary=factor(ti_binary,levels=c('yes','no'))), x='patient', y='active_RDI', fill='ti_binary', size=0.3, outlier.size=1) + 
  theme(legend.position = 'bottom', axis.title.x = element_blank(), axis.text.x = element_text(angle=-30,hjust=0))
fig5e = fig5e + stat_pvalue_manual(stat.test,  label = "p", tip.length = .02, y.position=9, bracket.size = 0.3, size=10/.pt)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels=label_number(accuracy=1))+
  scale_fill_manual(values=c('#D0759F','#88BECD'))+
  ylab('# of active RDIs')+
  labs(fill = "CD8+ T cell infiltration")
plot(fig5e)

## ------- Sup 9B Boxplot score between spatial spots with and without cd8 infiltration slide level -------
stat.test = res %>% mutate(ti=`TCD8 infiltration`) %>% group_by(name) %>% rstatix::wilcox_test(active_RDI ~ ti, ref.group = '1', alternative = 'greater')
stat.test <- stat.test %>%
  add_xy_position(x = "name", dodge = 0.8) %>% add_significance(p.col='p')
sup9b = ggboxplot(res%>% mutate(ti=`TCD8 infiltration`, ti_binary = ifelse(ti == 1,'yes','no')) %>% mutate(ti_binary=factor(ti_binary,levels=c('yes','no'))), x='name', y='active_RDI', fill='ti_binary', size=0.3, outlier.size=1) + 
  theme(legend.position = 'bottom', axis.title.x = element_blank(), axis.text.x = element_text(angle=-30,hjust=0))
sup9b = sup9b + stat_pvalue_manual(stat.test,  label = "p", tip.length = .02, y.position=9, bracket.size = 0.3, size=10/.pt)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels=label_number(accuracy=1))+
  scale_fill_manual(values=c('#D0759F','#88BECD'))+
  ylab('# of active RDIs')+
  labs(fill = "CD8+ T cell infiltration")
plot(sup9b)

## ------- Figure 5F Boxplot score between spatial spots in resistant vs naive spot level -------
stat.test = res %>% rstatix::wilcox_test(active_RDI ~ group, ref.group = 'naïve', alternative = 'greater')
stat.test <- stat.test %>%
  add_xy_position(x = "group", dodge = 0.8) %>% add_significance(p.col='p')

fig5f = ggviolin(res, x='group', y='active_RDI', fill='group', size=0.3) + 
  theme(legend.position = 'None', axis.title.x = element_blank(), axis.text.x = element_text(angle=-30,hjust=0))
fig5f = fig5f + stat_pvalue_manual(stat.test,  label = "p", tip.length = .02, y.position=9, bracket.size = 0.3, size=10/.pt)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels=label_number(accuracy=1))+
  scale_fill_manual(values=c('#D0759F','#88BECD'))+
  ylab('# of active RDIs')
plot(fig5f)

## ------- Figure 5G Boxplot score between spatial spots in resistant vs naive patient level -------
## slide specific may delete 
res_patientlvl = res %>% group_by(patient, group, name) %>% dplyr::summarise(active_RDI=sum(active_RDI)) %>% 
  ungroup() %>% group_by(patient, group) %>% dplyr::summarise(active_RDI=mean(active_RDI))
stat.test = wilcox.test(active_RDI ~ group, res_patientlvl, alternative='g')
stat.test = data.frame(group1=c('naïve'), group2=c('PD'), p=round(stat.test$p.value,2))

fig5g=ggplot(res_patientlvl, aes(x=group, y=active_RDI, fill=group, color=group), size=0.3)+
  geom_violin(trim=FALSE, fill=NA, show.legend = F)+
  geom_boxplot(width=0.3, color='black', show.legend = F)+
  theme_pubr()+
  xlab(NULL)+
  ylab('# of active RDIs')+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels=label_number(accuracy=1))+
  scale_fill_manual(values=c('#D0759F','#88BECD'))+
  scale_color_manual(values=c('#D0759F','#88BECD'))+
  theme(plot.title = element_text(hjust = 0.5),)+
  stat_compare_means(ref.group='naïve', paired=F, method='wilcox.test', method.args = list(alternative = "greater"), bracket.size = 0.4, label='p', 
                     comparisons = list(c('naïve','PD')), tip.length = 0.03, size=10/.pt, label.y=350)
plot(fig5g)

## ------- Figure 8AB Spatial Map -------

sup8a = ggplot(res, aes(x,y))+
  geom_point(aes(color=RDS_scale), size=0.5)+
  scale_color_viridis(option='magma', direction=-1)+
  facet_wrap(~name, nrow=4, ncol=4, strip.position = c('top'))+
  theme_pubr()+
  labs(color= '# of active RDIs (scaled)')
plot(sup8a)

sup8b=ggplot(res, aes(x,y))+
  geom_point(aes(color=`TCD8 fraction_scale`), size=.5)+
  scale_color_viridis(option='magma', direction=-1)+
  facet_wrap(~name, nrow=4, ncol=4)+
  theme_pubr()+
  labs(color= 'fraction of CD8+ T cells (scaled)')
plot(sup8b)

## ------ Calculate AUC to test overlap of pattern ------

df = res %>% mutate(ti = factor(`TCD8 infiltration`, levels=c(0,1)))

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
sd(auc_values)
median(auc_values)
