# Load required libraries

library(plyr)
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(ggplot2)
library(DataCombine)
library(parallel)
library(logistf)
library(Hmisc)
library(reshape2)
library(patchwork)
library(grid)
library(gridExtra)
library(gtable)
library(scales)
library(UpSetR)

args <- commandArgs(trailingOnly = TRUE)

source('Plot_functions.R')

start.time <- Sys.time()
print(start.time)

# Set required paths
scorefolder<-'results/construct_scores/'
analysisdir<- 'results/analysis_results/'
outcome_trained='senomi_phase4';outcome='senomi_phase4'
plots_dir='results/plots/'


if (!dir.exists(analysisdir)) {
  dir.create(analysisdir, recursive = TRUE)
}

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

geneticpredictors= c('clinicalvariant', 'gwastrait','geneburden','singlevar') 
geneticpredictorsdoe= data.frame(from = geneticpredictors, to = c("clinicalvariant_DOE","gwastrait_doe","geneburden_doe","singlevar_doe"))
Predictor_names= data.frame(predictors=geneticpredictors, New= c('ClinicalVariant','GWA trait', 'GeneBurden','SingleVar'))
Category_names= data.frame(predictors=c('ClinicalVariant','GWA trait', 'GeneBurden','SingleVar'), Category = c(rep('Clinical variants',1),rep('GWAS trait',1),rep('Coding variants',2)))


##############################################################################
# Figure 2 - Association of genetic features with drug side effects 
##############################################################################

# --- Run model ---

lapply(c('Opentargets','Onsides'), function(dataset){
  
  datafile<-list.files(path='Data/', pattern='dataset_se.mi_withgenetics', full.names=T)
  datafile=datafile[grepl(paste0(dataset), datafile)]
  dataset1<-fread(datafile, data.table=F)
  
  # Run univariate logistic regression models for each genetic predictor
  univarresults<-do.call(rbind,mclapply(geneticpredictors, function(pheno){
    
    # Identify all unique genes and parentterms where the predictor is non-zero and when the outcome is positive (is a side effect)
    genes_pred=dataset1 %>% distinct(gene,.data[[pheno]]) %>% filter(.data[[pheno]]!=0)
    PT_pred=dataset1 %>% distinct(parentterm,.data[[pheno]]) %>% filter(.data[[pheno]]!=0)
    genes_pred_outcome=dataset1 %>% distinct(gene,.data[[pheno]],.data[[outcome]]) %>% filter(.data[[pheno]]!=0) %>% filter(.data[[outcome]]==1)
    PT_pred_outcome=dataset1 %>% distinct(parentterm,.data[[pheno]],.data[[outcome]]) %>% filter(.data[[pheno]]!=0) %>% filter(.data[[outcome]]==1)
    
    # Fit logistic regression model with outcome ~ genetic predictor + category (adjusting for category)
    model1 =glm(as.formula(paste0(outcome, ' ~', pheno,'+category')), data=dataset1,family = 'binomial')

    # Extract regression results (coefficients, confidence intervals, p-values)
    ta1<-rbind(cbind.data.frame(outcome=outcome,Analysis='Full', predictors=names(model1$coefficients)[!grepl('category|Inter', names(model1$coefficients))],
                                No.genes=length(unique(genes_pred$gene)), No.parentterms=length(unique(PT_pred$parentterm)),
                                No.genes.outcome=length(unique(genes_pred_outcome$gene)), No.parentterms.outcome=length(unique(PT_pred_outcome$parentterm)),
                                OR=round(exp(summary(model1)$coefficient[,1])[!grepl('category|Inter', names(model1$coefficients))],3),
                                lowerCI=round(exp(summary(model1)$coefficient[,1][!grepl('category|Inter', names(model1$coefficients))]-(1.96* summary(model1)$coefficient[,2][!grepl('category|Inter', names(model1$coefficients))])),3),
                                upperCI=round(exp(summary(model1)$coefficient[,1][!grepl('category|Inter', names(model1$coefficients))]+(1.96* summary(model1)$coefficient[,2][!grepl('category|Inter', names(model1$coefficients))])),3),P.val=summary(model1)$coefficient[,4][!grepl('category|Inter', names(model1$coefficients))]))
    return(ta1)
  }, mc.cores=4))
  
  # --- Create forest plot ---
  
  options(digits = 3, scipen = -2)
  univar_data_outcome <- univarresults %>% mutate(upperCI=as.numeric(upperCI),Pval_sig=as.factor(ifelse(P.val>=0.05,0,1)))
  
  if(length(unique(univar_data_outcome$Pval_sig[!is.na(univar_data_outcome$Pval_sig)]))==1){
    shape_vals=c(19)
  } else{
    shape_vals=c(1,19)
  }
  #rename predictors for plot labels
  univar_data_outcome <- FindReplace(data = univar_data_outcome, Var = "predictors", replaceData = Predictor_names, from = "predictors", to = "New", exact = TRUE)
  univar_data_outcome <-inner_join(univar_data_outcome,Category_names) %>% arrange(Category)
  univar_data_outcome$order=rep(length(unique(univar_data_outcome$predictors)):1)
  univar_data_outcome$predictors<-factor(univar_data_outcome$predictors, levels=unique(univar_data_outcome$predictors[order(univar_data_outcome$order, decreasing=T)]))
  univar_data_outcome$CI= ifelse( univar_data_outcome$upperCI>100, paste0(round(univar_data_outcome$OR,1),' (', round(univar_data_outcome$lowerCI,1), ' - ', format(round(univar_data_outcome$upperCI,1), nsmall = 1),')')
                                  , paste0(round(univar_data_outcome$OR,1),' (', round(univar_data_outcome$lowerCI,1), ' - ', round(univar_data_outcome$upperCI,1),')'))
  ## Add OR as labels to plot
  ORlabels=univar_data_outcome %>% dplyr::mutate(Label=paste0( univar_data_outcome$CI))   %>% group_by(predictors) %>% dplyr::summarise(Label=paste0(Label, collapse='\n'))%>% select(predictors,Label)
  
  fp <- create_forest_plot(df = univar_data_outcome,
                           yvar = "predictors", 
                           groupvar = "Category",
                           shapevar = "Pval_sig",
                           colorvar = "Category",
                           ci_lower = "lowerCI",
                           ci_upper = "upperCI",
                           y_levels = univar_data_outcome$predictors,
                           label_above = paste0(round(univar_data_outcome$No.genes.outcome, 2),'/',round(univar_data_outcome$No.genes, 2)),
                           label_below = paste0(round(univar_data_outcome$No.parentterms.outcome ,3), '/',round(univar_data_outcome$No.parentterms ,3)),
                           shape_vals = shape_vals,
                           coord_limit = 4)
  
  p_right_all1 <- make_side_labels(df = ORlabels,
                           yvar = "predictors", 
                           label = "Label", 
                           y_levels = univar_data_outcome$predictors,
                           size = 3.1, 
                           hjust_size = 0.75)
  
  fp3 <- combine_plots(fp = fp,
                       labels_plot = p_right_all1, 
                       widths = c(4,-1.5,3.5),
                       legend = "bottom")
  
ggsave(fp3, file=paste0(plots_dir,'/Figure2_forestplot_univar_',dataset,'.png'), width = 6.5, height=5, dpi=300)
  
})


##############################################################################
# Figure 3 -Association of the SE-GPS with drug side effects by groupings
##############################################################################

# --- Prepare data: stratify dataset into n gene targets, n drugs per se and n se and by category ---
# loop through full dataset ('withdrugs') and severe drug dataset ('withdrawn') for both Open Targets and Onsides

lapply(c('withdrugs','withdrawn'), function(drugs){
  
  Stratified_Regression_analysis<-do.call(rbind, lapply(c('Onsides','Opentargets'), function(dataset){
    
    if(drugs=='withdrugs'){
      Genescores_beta_genetic_CV=fread(paste0(scorefolder, 'All_genescoresum_across_all_predictors_',dataset,'_outcome_',outcome,'_sideeffectproject.txt.gz'), data.table=F)
      Dataset_genescores1= Genescores_beta_genetic_CV %>%arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% mutate(percent=order/length(genescoresum) *100)
      datafile<-list.files(path='Data/', pattern='dataset_se.mi_withgenetics', full.names=T)
      datafile=datafile[grepl(paste0(dataset), datafile)]
      dataset1<-fread(datafile, data.table=F)
    } else{
      print(paste('withdrawn'))
      datafile<-list.files(path=paste0('Data/'), pattern='Withdrawndrugs_dataset_', full.names=T)
      datafile=datafile[grepl(paste0(dataset), datafile)]
      dataset1<-fread(datafile, data.table=F)
      Dataset_genescores1= dataset1 %>%arrange(genescoresum)
    }
    
    dataset1<-fread(datafile, data.table=F)
    
    # --- Summary Stats: Gene Targets, SE, Drugs ---
    
    # Number of gene targets per drug
    gene_targets <- Dataset_genescores1 %>%
      distinct(drugname, gene) %>% group_by(drugname) %>% tally() %>%
      rename(number_gene_targets = n) %>%
      mutate(number_gene_targets = case_when(number_gene_targets == 1 ~ '1',
                                             number_gene_targets %in% 2:5 ~ '2-5',
                                             number_gene_targets > 5 ~ '5+'))
    #how many SE per drug
    sefreq <- Dataset_genescores1 %>%
      distinct(drugname, se, parentterm) %>% filter(se==1) %>% group_by(drugname) %>% tally() %>%
      rename(number_SE = n) %>%
      mutate(number_SE = case_when(number_SE %in% 1:4 ~ '1-4',
                                   number_SE %in% 5:10 ~ '5-10',
                                   number_SE > 10 ~ '10+'))
    # how many drugs per SE ( drug specific )
    drugfreq <- Dataset_genescores1 %>%
      distinct(drugname, se, parentterm) %>% filter(se==1) %>% group_by(parentterm) %>% tally() %>%
      rename(number_drugs = n) %>%
      mutate(number_drugs = case_when(number_drugs %in% 1:10 ~ '1-10',
                                      number_drugs %in% 11:30 ~ '11-30',
                                      number_drugs > 30 ~ '30+'))
    
    Dataset_genescores2 <- Dataset_genescores1 %>%
      inner_join(gene_targets, by = "drugname") %>%
      inner_join(drugfreq, by = "parentterm") %>%
      inner_join(sefreq, by = "drugname") %>%
      mutate(
        number_SE = factor(number_SE, levels = c('1-4', '5-10', '10+')),
        number_gene_targets = factor(number_gene_targets, levels = c('1', '2-5', '5+')),
        number_drugs = factor(number_drugs, levels = c('1-10', '11-30', '30+')),
        score_all=ifelse(genescoresum>0,1,0)
      )
    
    # --- Run logistic regression with outcome as side effect and predictor as binarized SE-GPS  + category (adjusting for category) ---
    
    # --- Run mode across full dataset ---
    
    cov_pred=paste0('category')
    model1 =glm(as.formula(paste0(outcome, ' ~ score_all + ', cov_pred)), data=Dataset_genescores2,family = 'binomial')

        # Extract regression results (coefficients, confidence intervals, p-values)
    Fullmodel<-rbind(cbind.data.frame(Analysis='Full',
                                      No.drugs=length(unique(Dataset_genescores2$drugname)), No.genes=length(unique(Dataset_genescores2$gene)),
                                      No.parentterms=length(unique(Dataset_genescores2$parentterm)),
                                      No.outcome=length(unique(Dataset_genescores2$parentterm[Dataset_genescores2[outcome]==1])),
                                      No.genes.GE=length(unique(Dataset_genescores2$gene[Dataset_genescores2$genescoresum!=0])),
                                      No.parentterms.GE=length(unique(Dataset_genescores2$parentterm[Dataset_genescores2$genescoresum!=0])),
                                      outcome=outcome,outcome_trained=outcome_trained, dataset=dataset, Stratified_group=NA, Level=NA,category=NA,
                                      OR=exp(summary(model1)$coefficient[2,1]),
                                      lowerCI=exp(summary(model1)$coefficient[2,1]-(1.96* summary(model1)$coefficient[2,2])),
                                      upperCI=exp(summary(model1)$coefficient[2,1]+(1.96* summary(model1)$coefficient[2,2])),P.val=summary(model1)$coefficient[2,4]))
    
    # --- Run model for stratified data ---
    
    Stratified<-do.call(rbind,lapply(c('number_gene_targets','number_drugs','number_SE'), function(strat){
      gene_freq_gps<-do.call(rbind,lapply(sort(unique(Dataset_genescores2[[strat]])), function(strat_level){
       # Filter data for the stratification level 
        Dataset_genescores1_strat_gene2 <- Dataset_genescores2 %>% filter(.data[[strat]] == strat_level)
        # Filter to only SEs that occurred in stratified data
        Dataset_genescores1_stratspecific=Dataset_genescores1_strat_gene2 %>% filter(.data[[strat]]==strat_level) %>% filter(senomi_phase4==1)
        #Subset stratified data to those parentterms present in filtered side effect subset 
        Dataset_genescores1_strat_gene=subset(Dataset_genescores1_strat_gene2, (parentterm %in% Dataset_genescores1_stratspecific$parentterm)) 
        cov_pred=paste0('category')
        # Run logistic regression model and extract regression results
        model1 =glm(as.formula(paste0(outcome, ' ~ score_all + ', cov_pred)), data=Dataset_genescores1_strat_gene,family = 'binomial')
        ta1_gps_gene<-rbind(cbind.data.frame(Analysis='Stratified_group', No.drugs=length(unique(Dataset_genescores1_strat_gene$drugname)), No.genes=length(unique(Dataset_genescores1_strat_gene$gene)),
                                             No.parentterms=length(unique(Dataset_genescores1_strat_gene$parentterm)),
                                             No.outcome=length(unique(Dataset_genescores1_strat_gene$parentterm[Dataset_genescores1_strat_gene[outcome]==1])),
                                             No.genes.GE=length(unique(Dataset_genescores1_strat_gene$gene[Dataset_genescores1_strat_gene$genescoresum!=0])),
                                             No.parentterms.GE=length(unique(Dataset_genescores1_strat_gene$parentterm[Dataset_genescores1_strat_gene$genescoresum!=0])),
                                             outcome=outcome,outcome_trained=outcome_trained, dataset=dataset, Stratified_group=strat, Level=strat_level,category=NA,
                                             OR=exp(summary(model1)$coefficient[2,1]),
                                             lowerCI=exp(summary(model1)$coefficient[2,1]-(1.96* summary(model1)$coefficient[2,2])),
                                             upperCI=exp(summary(model1)$coefficient[2,1]+(1.96* summary(model1)$coefficient[2,2])),P.val=summary(model1)$coefficient[2,4]))
      }))
      return(gene_freq_gps)
    }))
    
    # --- Stratify by category ---
    
    Bycategory<-do.call(rbind,lapply(unique(Dataset_genescores2$category), function(categ){
      print(categ)
      # Filter data by phecode category
      Dataset_genescores1_cat1=Dataset_genescores2 %>% filter(categ==category)
      # Filter to only category SEs that occurred in stratified data
      Dataset_genescores1_catspec=Dataset_genescores2 %>% filter(categ==category) %>% filter(senomi_phase4==1)
      # Subset to drugs with SE that match phecode category
      Dataset_genescores1_cat=subset(Dataset_genescores1_cat1, (drugname %in% Dataset_genescores1_catspec$drugname)) 
      # Make sure score and side effect have two values to avoid model failing
      if(length(unique(Dataset_genescores1_cat$score_all))>1 & length(unique(Dataset_genescores1_cat$senomi_phase4))>1) {
        model1 =glm(as.formula(paste0(outcome, ' ~ score_all ')), data=Dataset_genescores1_cat,family = 'binomial')
        ta1_gps_cat<-rbind(cbind.data.frame( Analysis='ByCategory',
                                             No.drugs=length(unique(Dataset_genescores1_cat$drugname)), No.genes=length(unique(Dataset_genescores1_cat$gene)),
                                             No.parentterms=length(unique(Dataset_genescores1_cat$parentterm)),
                                             No.outcome=length(unique(Dataset_genescores1_cat$parentterm[Dataset_genescores1_cat[outcome]==1])),
                                             No.genes.GE=length(unique(Dataset_genescores1_cat$gene[Dataset_genescores1_cat$genescoresum!=0])),
                                             No.parentterms.GE=length(unique(Dataset_genescores1_cat$parentterm[Dataset_genescores1_cat$genescoresum!=0])),
                                             
                                             outcome=outcome,outcome_trained=outcome_trained,dataset=dataset,Stratified_group=NA, Level=NA, category=categ,
                                             OR=exp(summary(model1)$coefficient[2,1]),
                                             lowerCI=exp(summary(model1)$coefficient[2,1]-(1.96* summary(model1)$coefficient[2,2])),
                                             upperCI=exp(summary(model1)$coefficient[2,1]+(1.96* summary(model1)$coefficient[2,2])),P.val=summary(model1)$coefficient[2,4]))
      }
    }))
    Stratified_analysis<-rbind(Fullmodel, Stratified,Bycategory)
    
  }))
  write.table(Stratified_Regression_analysis, paste0(analysisdir,'Stratified_analysis_full_outcome_trained',outcome_trained,'_',drugs,'.txt'), sep='\t', quote=F, row.names=F)
  
})

# --- Create forest plot ---

lapply(c('withdrugs','withdrawn'), function(drugs){
  
  lapply(c('Opentargets','Onsides'), function(dataset){
    
    outcomestrat<-fread(paste0(analysisdir,'Stratified_analysis_full_outcome_trained',outcome_trained,'_',drugs,'.txt'), data.table=F)
    univar_data_outcome1 <- outcomestrat %>% filter( dataset == !!dataset)
    # Format plot labels
    univar_data_outcome1 <- univar_data_outcome1 %>% dplyr::mutate(Stratified_group=gsub('_','',gsub('number_','N ',Stratified_group )), Stratified_group= gsub('genetargets','gene targets',Stratified_group),Stratified_group= gsub('N drugs','N drugs per SE',Stratified_group))
    univar_data_outcome1$Predictor=paste0(univar_data_outcome1$Stratified_group, ': ', univar_data_outcome1$Level, sep = ' ')
    univar_data_outcome1$Predictor=gsub('NA: NA ','All',univar_data_outcome1$Predictor)
    univar_data_outcome1$grouppred=gsub(':.*','',univar_data_outcome1$Predictor)
    univar_data_outcome_drugs = univar_data_outcome1 %>% filter(grouppred %in% c('All','N gene targets','N drugs per SE','N SE' )) %>% filter(is.na(category))
    univar_data_outcome_drugs$Pval_sig=ifelse(univar_data_outcome_drugs$P.val>=0.05,0,1)
    univar_data_outcome_drugs$Pval_sig<-as.factor(univar_data_outcome_drugs$Pval_sig)
    
    if(length(unique(univar_data_outcome_drugs$Pval_sig[!is.na(univar_data_outcome_drugs$Pval_sig)]))==1){
      shape_vals=c(19)
    } else{
      shape_vals=c(1,19)
    }
    univar_data_outcome_drugs$order=rep(length(unique(univar_data_outcome_drugs$Predictor)):1)
    
    univar_data_outcome_drugs$Predictor<-factor(univar_data_outcome_drugs$Predictor, levels=unique(univar_data_outcome_drugs$Predictor))
    maxor=ifelse(max(univar_data_outcome_drugs$upperCI[!is.na(univar_data_outcome_drugs$upperCI)])>10, 10, max(univar_data_outcome_drugs$upperCI[!is.na(univar_data_outcome_drugs$upperCI)]))
    univar_data_outcome_drugs$upperCI_cut<-ifelse(univar_data_outcome_drugs$upperCI> maxor, maxor, NA)
    univar_data_outcome_drugs$CI= paste0(round(univar_data_outcome_drugs$OR,1),' (', round(univar_data_outcome_drugs$lowerCI,1), ' - ', round(as.numeric(univar_data_outcome_drugs$upperCI),1),')')
    
    ORlabels=univar_data_outcome_drugs %>% mutate(Label=paste0(univar_data_outcome_drugs$CI)) %>% group_by(Predictor)  %>% select(Predictor,Label)
    
    fp <- create_forest_plot(df = univar_data_outcome_drugs,
                             yvar = "Predictor", 
                             shapevar = "Pval_sig",
                             colorvar = "grouppred",
                             ci_lower = "lowerCI",
                             ci_upper = "upperCI",
                             y_levels = univar_data_outcome_drugs$Predictor,
                             label_above = paste(round(univar_data_outcome_drugs$No.genes, 2)),
                             label_below = paste(round(univar_data_outcome_drugs$No.parentterms ,3)),
                             shape_vals = shape_vals,
                             coord_limit = max(univar_data_outcome_drugs$upperCI))
    
    # if CI are too long add arrows at upper cut off 
    if (any(!is.na(univar_data_outcome_drugs$upperCI_cut))) {
      fp <- add_arrows(fp, maxor, "Predictor")
    }
    
    p_right_all <- make_side_labels(df = ORlabels,
                                     yvar = "Predictor", 
                                     label = "Label", 
                                     y_levels = univar_data_outcome_drugs$Predictor,
                                     size = 3.1, 
                                     hjust_size = 0.72 )
    
    fp3 <- combine_plots(fp = fp,
                         labels_plot = p_right_all, 
                         widths = c(4,-1.5,3.5),
                         legend = "none")
    
    
    if(drugs=='withdrugs'){
      filename1=paste0(plots_dir,'/Figure3a_forestplot_stratified_GPS_',dataset,'.png')
    } else {
      filename1=paste0(plots_dir,'/Sup_Fig10a_forestplot_stratified_GPS_',dataset,'.png')
    }
    
    ggsave(fp3, file=filename1, width = 6.5, height=5)
    
    # --- Create forest plot for categories ---
    
    outcomestrat<-fread(paste0(analysisdir,'Stratified_analysis_full_outcome_trained',outcome_trained,'_',drugs,'.txt'), data.table=F)
    outcomestrat_bycat<-outcomestrat %>% filter(Analysis=='ByCategory'|Analysis=='Full', dataset == !!dataset)
    phecodef=read.csv('Data/phecodeX_info.csv')
    phecodef = phecodef %>% dplyr::rename(category_full=category) %>% mutate(category=gsub('_.*', '',phecode)) %>% select(category_full, category ) %>% distinct()
    outcomestrat_bycat1=left_join(outcomestrat_bycat,phecodef )
    outcomestrat_bycat1$category_full[outcomestrat_bycat1$Analysis=='Full']<-'All'
    outcomestrat_bycat1$upperCI<- as.numeric(outcomestrat_bycat1$upperCI)
    outcomestrat_bycat1$Pval_sig<-as.factor(ifelse(outcomestrat_bycat1$P.val>=0.05,0,1))
    
    if(length(unique(outcomestrat_bycat1$Pval_sig[!is.na(outcomestrat_bycat1$Pval_sig)]))==1){
      shape_vals=c(19)
    } else{
      shape_vals=c(1,19)
    }
    
    maxor=max(outcomestrat_bycat1$OR)+2
    outcomestrat_bycat1$upperCI_cut<-ifelse(outcomestrat_bycat1$upperCI> maxor+2.5 & outcomestrat_bycat1$upperCI!=Inf, maxor+2.5, NA)
    outcomestrat_bycat1$CI= paste0(round(outcomestrat_bycat1$OR,1),' (', round(outcomestrat_bycat1$lowerCI,1), ' - ', round(as.numeric(outcomestrat_bycat1$upperCI),1),')')
    outcomestrat_bycat1$upperCI[outcomestrat_bycat1$upperCI==Inf]<-0
    outcomestrat_bycat1=outcomestrat_bycat1 %>% arrange(desc(Analysis),desc(OR))
    outcomestrat_bycat1$order=rep(length(unique(outcomestrat_bycat1$category_full)):1)
    outcomestrat_bycat1$category_full<-factor(outcomestrat_bycat1$category_full, levels=unique(outcomestrat_bycat1$category_full))
    outcomestrat_bycat1$upperCI=ifelse(!is.na(outcomestrat_bycat1$upperCI_cut),outcomestrat_bycat1$upperCI_cut,outcomestrat_bycat1$upperCI)
    
    ORlabels=outcomestrat_bycat1 %>% mutate(Label=paste0(CI)) %>% select(category_full,dataset,Label) %>% arrange(desc(dataset)) %>% group_by(category_full) %>% dplyr::summarise(Label=paste0(Label,collapse='\n') )
    
    if(drugs=='withdrawn'){
      maxor = maxor + 4
    }
    fp <- create_forest_plot(df = outcomestrat_bycat1,
                             yvar = "category_full", 
                             shapevar = "Pval_sig",
                             colorvar = "category_full",
                             ci_lower = "lowerCI",
                             ci_upper = "upperCI",
                             y_levels = outcomestrat_bycat1$category_full,
                             shape_vals = shape_vals)
    
    # if CI are too long add arrows at upper cut off 
    if (any(!is.na(outcomestrat_bycat1$upperCI_cut))) {
      fp <- add_arrows(fp, maxor, "category_full")
    }
    
    p_right_all <- make_side_labels(df = ORlabels,
                                    yvar = "category_full", 
                                    label = "Label", 
                                    y_levels = outcomestrat_bycat1$category_full,
                                    size = 3.1, 
                                    hjust_size = 0.72 )
    
    fp3 <- combine_plots(fp = fp,
                         labels_plot = p_right_all, 
                         widths = c(4,-1.5,3.5),
                         legend = "none")
    

    if(drugs=='withdrugs'){
      filename1=paste0(plots_dir,'/Figure3c_forestplot_stratifiedbycategory_GPS_',dataset,'.png')
    } else {
      filename1=paste0(plots_dir,'/Sup_Fig10c_forestplot_stratifiedbycategory_GPS_',dataset,'.png')
    }
    ggsave(fp3, file=filename1, width = 7, height=6)
    
  })
})


##############################################################################
# Figure 4 and 5 - Association of the SE-GPS with drug side effects 0.3 bins
##############################################################################

# --- Prepare data and run model: loop across 0.3 increment bins  --- 
returnunivar<-do.call(rbind,lapply(c('Onsides','Opentargets'), function(dataset){
  do.call(rbind,lapply(c('withdrugs','withdrawn'), function(drugs){
    if(drugs=='withdrugs'){
      Genescores_beta_genetic_CV=fread(paste0(scorefolder, 'All_genescoresum_across_all_predictors_',dataset,'_outcome_',outcome,'_sideeffectproject.txt.gz'), data.table=F)
      Dataset_genescores1= Genescores_beta_genetic_CV %>%arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% mutate(percent=order/length(genescoresum) *100)
    } else{
      print(paste('withdrawn'))
      datafile<-list.files(path=paste0('Data/'), pattern='Withdrawndrugs_dataset_', full.names=T)
      datafile=datafile[grepl(paste0(dataset), datafile)]
      dataset1<-fread(datafile, data.table=F)
      Dataset_genescores1= dataset1 %>%arrange(genescoresum)
    }
    
    cov_pred=paste0('category')
    
    # calculate the number of non-zero genetic predictors per row and store in new column 'no.pred'
    Dataset_genescores1$no.pred=rowSums(Dataset_genescores1[,geneticpredictors] != 0)
    Dataset_genescores1[,geneticpredictors][Dataset_genescores1[,geneticpredictors] !=0 ] <- 1
    ##  Loop across thresholds of SE-GPS (from 0 to 2.1 by increments of 0.3)
    binned_gps2<-do.call(rbind,mclapply(c(seq(0,2.1,0.3)), function(bin) {
      # Filter dataset for entries where the SE-GPS is greater than the current threshold and binarize the predictor as 'bin'
      df_filtered=Dataset_genescores1 %>% filter(genescoresum>(bin) ) %>% mutate(bin=1)
      # Filter dataset with no genetic evidence
      Dataset_genescores1_below=Dataset_genescores1 %>% filter(no.pred==0 ) %>% mutate(bin=0)
      Dataset_genescorescompare=rbind(df_filtered , Dataset_genescores1_below) %>% distinct()
      # Get counts of genes, phenotypes with score > bin increment 
      genes=length(unique(Dataset_genescorescompare$gene[Dataset_genescorescompare$bin==1]))
      parentterm=length(unique(Dataset_genescorescompare$parentterm[Dataset_genescorescompare$bin==1]))
      gene_pheno=nrow(df_filtered %>% distinct(gene, parentterm))/nrow(Dataset_genescores1 %>% distinct(gene, parentterm)) *100
      entries=nrow(Dataset_genescorescompare[Dataset_genescorescompare$bin==1,])
      # Fit logistic regression model with increment threshold bin as predictor and include category as a covariate
      model1 =glm(as.formula(paste0(outcome,' ~ bin + category')), data=Dataset_genescorescompare,family = 'binomial')
      mod_output<-rbind(cbind.data.frame(dataset=dataset, outcome=outcome,binseq=0.3, drugs=drugs,genescoresum =bin, Percentile=round(df_filtered$percent[1],4) , genes=genes, parentterms=parentterm,gene_pheno_percentage=gene_pheno, entries=entries, OR=exp(summary(model1)$coefficient[2,1]),lowerCI=exp(summary(model1)$coefficient[2,1]-(1.96* summary(model1)$coefficient[2,2])),upperCI=exp(summary(model1)$coefficient[2,1]+(1.96* summary(model1)$coefficient[2,2])),P.val=summary(model1)$coefficient[2,4], analysis=paste0('All')))
      return(mod_output)
    }, mc.cores=10))
    binned_gps2$CI=paste0(round(binned_gps2$lowerCI,1), '-',round(binned_gps2$upperCI,1))
    binned_gps2$P.val=ifelse(binned_gps2$P.val==0,round(.Machine$double.xmin,311),binned_gps2$P.val)
    
    write.table(binned_gps2, paste0(analysisdir,'Binned_by_03_increments_',drugs,'_', dataset,'.txt'), sep='\t', quote=F, row.names=F)
    
  }))
}))

# --- Create increased bin plot ---

lapply(c('Onsides','Opentargets'), function(dataset){
  lapply(c('withdrugs','withdrawn'), function(drugs){
    
    binned_gps2=fread(paste0(analysisdir,'Binned_by_03_increments_',drugs,'_', dataset,'.txt'), data.table=F)
    binned_gps2$upperCI=as.numeric(binned_gps2$upperCI)
    max_value <- round(max(binned_gps2$OR)) + 5
    binned_gps2$upperCI_cut<-ifelse(binned_gps2$upperCI>max_value ,max_value,NA)
    binned_gps2$or_sig<-as.factor(ifelse(binned_gps2$P.val>=0.05,0,1))
    binned_gps2$genescoresum=factor(binned_gps2$genescoresum)
    if(length(unique(binned_gps2$or_sig[!is.na(binned_gps2$or_sig)]))==1){
      shape_vals=c(19)
    } else{
      shape_vals=c(1,19)
    }
    
    max_uci <- round(ifelse(max(binned_gps2$upperCI[!is.na(binned_gps2$upperCI)])>max_value,max_value , max(binned_gps2$upperCI[!is.na(binned_gps2$upperCI)])+1))
    
    plot_binned_or_with_table(
      binned_gps2 = binned_gps2,
      dataset = dataset, 
      drugs = drugs,
      title_PLT = "",
      score_plt = "score1",
      shape_vals = shape_vals,
      max_uci = max_uci,
      filename = if(drugs=='withdrugs'){
        filename1=paste0(plots_dir,'/Figure4_Increase_GPS_binned03_',dataset,'_',drugs,'_withtable','_prac.png')
      } else {
        filename1=paste0(plots_dir,'/Figure5_Increase_GPS_binned03_',dataset,'_',drugs,'_withtable','_prac.png')
      }
    )
    
    
  })
})

##############################################################################
# Figure 6 -Association of the SE-GPS-DOE with drug side effects 
##############################################################################

UNIVAR<-do.call(rbind,lapply(c('Onsides','Opentargets'), function(dataset){
  #Load SE-GPS-DOE data
  Genescores_beta_genetic_CV=fread(paste0(scorefolder, 'All_genescoresum_across_all_predictors_',dataset,'_outcome_',outcome_trained,'_sideeffectproject_DOE.txt.gz'), data.table=F)
  doedatafile<-list.files(path='Data/', pattern='DOE', full.names=T)
  doedatafile=doedatafile[grepl(paste0(dataset), doedatafile)]
  print(doedatafile)
  # Annotate scores as 'pos', 'neg', or 0
  Genescores_beta_genetic_CV$score=0
  Genescores_beta_genetic_CV$score= ifelse(Genescores_beta_genetic_CV$genescoresum >0,  'pos' , Genescores_beta_genetic_CV$score)
  Genescores_beta_genetic_CV$score= ifelse(Genescores_beta_genetic_CV$genescoresum <0,  'neg' , Genescores_beta_genetic_CV$score)
  dataset1<-fread(doedatafile, data.table=F)
  Dataset_genescores1= Genescores_beta_genetic_CV %>%arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% mutate(percent=order/length(genescoresum) *100)
   # Loop over positive and negative scores separately
  DOE_assoc<-do.call(rbind,lapply(c('pos','neg'), function(score){
    if(score=='pos') {
      Genescores_beta_genetic_CV$genescoresum_pos= ifelse(Genescores_beta_genetic_CV$genescoresum >0,  Genescores_beta_genetic_CV$genescoresum , 0)
      Dataset_genescores1= Genescores_beta_genetic_CV%>% select(-genescoresum) %>% dplyr::rename(genescoresum=genescoresum_pos) %>% arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% mutate(percent=order/length(genescoresum) *100)
    }else {
      Genescores_beta_genetic_CV$genescoresum_neg= ifelse(Genescores_beta_genetic_CV$genescoresum < 0,  Genescores_beta_genetic_CV$genescoresum , 0)
      Dataset_genescores1= Genescores_beta_genetic_CV%>% select(-genescoresum) %>% dplyr::rename(genescoresum=genescoresum_neg) %>% mutate(genescoresum=genescoresum*-1) %>% arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% mutate(percent=order/length(genescoresum) *100)
    }
    # create binary outcome: 1 if SE-GPS-DOE > 0
    Dataset_genescores1$score_all=ifelse(Dataset_genescores1$genescoresum>0,1,0)
    cov_pred=paste0('category') 
    # Set the MOA (mechanism of action) direction
    if(score=='pos') {moadrug='inhibitor'}
    if(score=='neg') {moadrug='activator'}
    
    # restrict to either activator or inhbitor drugs and parentterms that have a side effect
    dataset1_moa=dataset1 %>% filter(moa==moadrug)
    dataset1_moa_pt=dataset1_moa %>% filter(se==1)
    dataset1_moa1=subset(dataset1_moa, (parentterm %in% dataset1_moa_pt$parentterm))
    dataset1_moa1= dataset1_moa1 %>% select(drugname:moa)
    Dataset_genescores_score=Dataset_genescores1 %>% select(drugname , gene, parentterm,score_all)
    # join with SE-GPS-DOE data
    Dataset_genescores2=inner_join(Dataset_genescores_score,dataset1_moa1)
    
    cov_pred=paste0('category') 
    # Fit logistic regression with side effect as outcome and binarized score as predictor with category as covariate
    model1 =glm(as.formula(paste0(outcome, ' ~ score_all + ', cov_pred)), data=Dataset_genescores2,family = 'binomial')
    ta1_gps_direction<-rbind(cbind.data.frame( No.outcome=length(unique(Dataset_genescores2$parentterm[Dataset_genescores2[outcome]==1])),
                                               outcome=outcome,dataset=dataset,score=score,moa=moadrug, 
                                               OR=exp(summary(model1)$coefficient[2,1]), 
                                               lowerCI=exp(summary(model1)$coefficient[2,1]-(1.96* summary(model1)$coefficient[2,2])),
                                               upperCI=exp(summary(model1)$coefficient[2,1]+(1.96* summary(model1)$coefficient[2,2])),P.val=summary(model1)$coefficient[2,4]))
    
    
    ## binned analysis 
    cols=c("clinicalvariant_DOE","gwastrait_doe","geneburden_doe","singlevar_doe")
    ## calculate number of genetic predictors present 
    Dataset_genescores1$no.pred=rowSums(Dataset_genescores1[,cols] != 0)  
    Dataset_genescores1[,cols][Dataset_genescores1[,cols] !=0 ] <- 1
   
     # Prepare data for binned analysis (separately for pos scores with inhibitor moa and neg scores with activator moa)
      dataset1_moa=dataset1 %>% filter(moa==moadrug)
      dataset1_moa_pt=dataset1_moa %>% filter(se==1)
      dataset1_moa1=subset(dataset1_moa, (parentterm %in% dataset1_moa_pt$parentterm))
      dataset1_moa1= dataset1_moa1 %>% select(drugname:moa)
      Dataset_genescores_score=Dataset_genescores1 %>% select(drugname , gene, parentterm,genescoresum,no.pred,percent)
      
      Dataset_genescores2=inner_join(Dataset_genescores_score,dataset1_moa1)
      cov_pred=paste0('category') 
      
      ##  Loop across thresholds of SE-GPS-DOE (from 0 to 2.1 by increments of 0.3)
      binned_gps3<-do.call(rbind,mclapply(c(seq(0,2.1,0.3)), function(bin) {
        print(bin)
        # Filter dataset for entries where the SE-GPS-DOE is greater than the current threshold and binarize the predictor as 'bin'
        df_filtered=Dataset_genescores2 %>% filter(genescoresum>=(bin) ) %>% mutate(bin=1)
        Dataset_genescores1_below=Dataset_genescores2 %>% filter(no.pred==0 ) %>% mutate(bin=0)
        Dataset_genescores6=rbind(df_filtered , Dataset_genescores1_below) %>% distinct()
        # Summary
        genes=length(unique(Dataset_genescores6$gene[Dataset_genescores6$bin==1]))
        parentterm=length(unique(Dataset_genescores6$parentterm[Dataset_genescores6$bin==1]))
        gene_pheno=nrow(df_filtered %>% distinct(gene, parentterm))/nrow(Dataset_genescores2 %>% distinct(gene, parentterm)) *100
        entries=nrow(Dataset_genescores6[Dataset_genescores6$bin==1,])
        # Logistic regression per bin
        model1 =glm(as.formula(paste0(outcome,' ~ bin +',cov_pred )), data=Dataset_genescores6,family = 'binomial')
        mod_output<-rbind(cbind.data.frame(dataset=dataset, outcome=outcome,score=score, genescoresum =bin, Percentile=round(df_filtered$percent[1],5) , genes=genes, parentterms=parentterm,gene_pheno_percentage=gene_pheno, entries=entries, OR=exp(summary(model1)$coefficient[2,1]),lowerCI=exp(summary(model1)$coefficient[2,1]-(1.96* summary(model1)$coefficient[2,2])),upperCI=exp(summary(model1)$coefficient[2,1]+(1.96* summary(model1)$coefficient[2,2])),P.val=summary(model1)$coefficient[2,4],
                                           analysis=paste0(moadrug)))
        return(mod_output)   
      }, mc.cores=10))                      
      
      binned_gps3 %>% select(genescoresum, Percentile, OR,entries) %>% mutate(percent=100-Percentile, OR=round(OR,2))
      binned_gps3$CI=paste0(round(binned_gps3$lowerCI,1), '-',round(binned_gps3$upperCI,1))
      binned_gps3$P.val=ifelse(binned_gps3$P.val==0,round(.Machine$double.xmin,311),binned_gps3$P.val)
      
    write.table(binned_gps3, paste0(analysisdir,'Binned_by_sum_binsize0.3_',dataset,'_',score,'_DOE_',moadrug,'.txt'), sep='\t', quote=F, row.names=F)
    
    return(ta1_gps_direction)
  }))
}))
# Save univariate analysis results

write.table(UNIVAR, paste0(analysisdir,'Univar_segps_doe_DOE_moadrug.txt'), sep='\t', quote=F, row.names=F)

# --- Create plot ---

UNIVAR<-fread(paste0(analysisdir,'Univar_segps_doe_DOE_moadrug.txt'), data.table=F)
# format data for plot
univar_doe <- UNIVAR %>%
  dplyr::mutate(
    score = gsub('pos', 'Positive\nSE-GPS DOE', score),
    score = gsub('neg', 'Negative\nSE-GPS DOE', score),
    Pval_sig = as.factor(ifelse(P.val >= 0.05, 0, 1)),
    order = as.numeric(factor(score, levels = rev(unique(score))))
  )

if(length(unique(univar_doe$Pval_sig[!is.na(univar_doe$Pval_sig)]))==1){
  shape_vals=c(19)
} else{
  shape_vals=c(1,19)
}
univar_doe$score<-factor(univar_doe$score, levels=unique(univar_doe$score))
maxor=ifelse(max(univar_doe$upperCI[!is.na(univar_doe$upperCI)])>8, 8, max(univar_doe$upperCI[!is.na(univar_doe$upperCI)]))
univar_doe$upperCI_cut<-ifelse(univar_doe$upperCI> maxor, maxor, NA)
univar_doe$CI= paste0(round(univar_doe$OR,1),' (', round(univar_doe$lowerCI,1), ' - ', round(as.numeric(univar_doe$upperCI),1),')')
univar_doe$dataset=ifelse(univar_doe$dataset=='Onsides','Onsides', 'Open Targets')
univar_doe$dataset <- factor(univar_doe$dataset, levels = (unique(univar_doe$dataset)))

ORlabels=univar_doe  %>% select(score,dataset,CI) %>% arrange(desc(dataset)) %>% group_by(score) %>% dplyr::summarise(Label=paste0(CI,collapse='\n') )

fp <- create_forest_plot(df = univar_doe,
                         yvar = "score", 
                         shapevar = "Pval_sig",
                         colorvar = "dataset",
                         ci_lower = "lowerCI",
                         ci_upper = "upperCI",
                         y_levels = univar_doe$score,
                         shape_vals = shape_vals,
                         coord_limit = max(univar_doe$upperCI))

  # if CI are too long add arrows at upper cut off 
if (any(!is.na(univar_doe$upperCI_cut))) {
    fp <- add_arrows(fp, maxor, "score")
  }

                         

p_right_all <- make_side_labels(df = ORlabels,
                                yvar = "score", 
                                label = "Label", 
                                y_levels = univar_doe$score,
                                size = 3.1, 
                                hjust_size = 0.72 )

fp3 <- combine_plots(fp = fp,
                     labels_plot = p_right_all, 
                     widths = c(4,-1.5,3.5),
                     legend = "bottom")

Filename1=paste0(plots_dir,'/Figure6_univar_DOE_moamatch.png')
ggsave(fp3, file=Filename1, width = 5, height=4)

# --- Create binned by increment se-gps-doe plot ---

lapply(c('Onsides','Opentargets'), function(dataset){
  
  lapply(c('pos','neg'), function(score_plt){
    if(score_plt=='pos'){
      binned_gps1=fread(paste0(analysisdir,'Binned_by_sum_binsize0.3_',dataset,'_',score_plt,'_DOE_inhibitor.txt'), data.table=F)
    } else{
      binned_gps1=fread(paste0(analysisdir,'Binned_by_sum_binsize0.3_',dataset,'_',score_plt,'_DOE_activator.txt'), data.table=F)
    }
    binned_gps2=binned_gps1[!is.na(binned_gps1$Percentile),]
    binned_gps2=binned_gps2  %>% filter(genescoresum!=2.1)  
    max_value=20
    binned_gps2$upperCI_cut<-ifelse(binned_gps2$upperCI>max_value ,binned_gps2$OR+6,NA) 
    binned_gps2$or_sig<-ifelse(binned_gps2$P.val>=0.05,0,1)
    binned_gps2$or_sig<-as.factor(binned_gps2$or_sig)
    binned_gps2$genescoresum=factor(binned_gps2$genescoresum)
    
    if(length(unique(binned_gps2$or_sig[!is.na(binned_gps2$or_sig)]))==1){
      shape_vals=c(19)
    } else{
      shape_vals=c(1,19)
    }
    rep_str = c('abs'='Side effect score\nDOE','neg'='Negative SE-GPS-DOE','pos'='Positive SE-GPS-DOE')
    binned_gps2$DOE <- str_replace_all(binned_gps2$score, rep_str)
    binned_gps2$DOE<-factor(binned_gps2$DOE, levels=unique(binned_gps2$DOE))
    
    binned_gps2$upperCI=as.numeric(binned_gps2$upperCI)
    max_uci <- round(ifelse(max(binned_gps2$upperCI[!is.na(binned_gps2$upperCI)])>max_value,max_value , max(binned_gps2$upperCI[!is.na(binned_gps2$upperCI)])+1))
    
    title_PLT=unique(binned_gps2$DOE)
    
    plot_binned_or_with_table(
      binned_gps2 = binned_gps2,
      dataset = dataset, 
      title_PLT = title_PLT, 
      shape_vals = shape_vals,
      max_uci = max_uci,
      filename=if(dataset=='Opentargets'){
        paste0(plots_dir,'/Sup_Fig16_Increase_SEGPS_DOE_binned03_',dataset,'_',score_plt,'_prac.png')
      } else{
        paste0(plots_dir,'/Sup_Fig17_Increase_SEGPS_DOE_binned03_',dataset,'_',score_plt,'_prac.png')
      }
        )
  })
})


#############################################
# Create additional Supplementary figures
#############################################


#  Supplementary Fig. 1  and Supplementary Fig. 2

library(UpSetR)
library(viridis)

lapply(c('All','lof','gof'),function(analysis){
  if(analysis=='All'){
    phenotype_allgene_pt2<-fread(paste0('Data/Allgenes_dataset_parentterm_predictors_collapsed_sideeffect_project.txt.gz'), data.table=F)
    phenotype_allgene_pt2= phenotype_allgene_pt2 %>% filter(!category_full %in% c('Neonatal','Pregnancy'))
    
    #### All predictors
    rep_str = c('omim'='OMIM','hgmd'='HGMD','clinvar'='ClinVar','ravargb'='RAVAR_GB','OTgb'='OT_GB',
                'ravarsv'='RAVAR_SV','genebasssv'= 'Genebass_SV','l2g'='L2G','eqtlpheno'='eQTL phenotype')
    colnames(phenotype_allgene_pt2) <- str_replace_all(colnames(phenotype_allgene_pt2), rep_str)
    
    df_pt1<-phenotype_allgene_pt2 %>% select(-category_full, -category,-gene,-parentterm, ) 
    filename=paste0(plots_dir,'/Sup_Fig2_descriptiveplt_upset_plot_allgenes.png')
    
  } else{
    rep_str = c('omim_predicted'='OMIM','hgmd_predicted'='HGMD','clinvar_predicted'='ClinVar','ravar_genburden_DOE'='RAVAR_GB','OTgeneburden_doe'='OT_GB',
                'ravar_snp_predicted'='RAVAR_SV','genebass_predicted'= 'Genebass_SV','OTG_DOE'='L2G','eqtl_Pheno_DOE'='eQTL phenotype')
    DOE=fread(paste0('Data/Predictors_with_direction_effect_sideeffect_project_',analysis,'_predvalue.txt'),data.table=F)
    DOE$parentterm=as.character(trunc(as.numeric(DOE$Phecode)))
    DOE1 <- aggregate(
      cbind(omim_predicted, hgmd_predicted, clinvar_predicted, genebass_predicted, ravar_snp_predicted,
            OTgeneburden_doe, ravar_genburden_DOE,OTG_DOE,eqtl_Pheno_DOE ) ~ gene + parentterm + category_full, data = DOE,FUN = max,na.rm = TRUE)
    
    DOE1[, 4:12] <- lapply(DOE1[, 4:12], function(x) ifelse(x !=0, 1, x))
    colnames(DOE1) <- str_replace_all(colnames(DOE1), rep_str)
    
    df_pt1<-DOE1 %>% select(-category_full,-gene,-parentterm)
    df_pt1[is.na(df_pt1)] <- 0
    df_pt1[] <- lapply(df_pt1[], function(x) as.numeric(as.character(x)))
    if(analysis=='lof') { filename=paste0(plots_dir,'/Sup_Fig11_descriptiveplt_upset_plot_lofcounts.png')}
    if(analysis=='gof') { filename=paste0(plots_dir,'/Sup_Fig12_descriptiveplt_upset_plot_gofcounts.png')}
    
  }
  
  png(filename, res=300, width =8, height=6, units='in' )
  print(upset(df_pt1, 
              nintersects = 25, 
              nsets = length(names(df_pt1)), 
              order.by = "freq", 
              decreasing = T, 
              mb.ratio = c(0.6, 0.4),
              number.angles = 20, 
              text.scale = c(1.3, 1.1, 1.3, 1.1, 1.3, 0.8), 
              point.size = 2.8, 
              line.size = 1.7,
              matrix.color="#1f78b4",
              main.bar.color="#1f78b4", 
              sets.bar.color="#1f78b4",
              mainbar.y.label = "Number of gene-phenotype entries",
              sets.x.label = "Total observations"
  ))
  dev.off()
})

# Supplementary Fig 3

df_counts_percat<-as.data.frame(do.call(rbind,lapply(c('Opentargets','Onsides'), function(dataset){
  df=fread(paste0('Data/',dataset,'_dataset_se.mi_withgeneticsfiltered_5_All_parentterm_drugse_final.txt.gz'), data.table=F)
  SE=df %>% filter(senomi_phase4==1) %>%distinct(drugname, parentterm,category_full) 
  SE1= SE %>% group_by(category_full) %>% tally()
  total_count <- sum(SE1$n)
  SE1 <- SE1 %>%
    mutate(percentage = (n / total_count) * 100) %>% mutate(dataset=dataset)
  return(SE1)
})))
df_counts_percat$order=rep(length(unique(df_counts_percat$category_full)):1)
df_counts_percat$category_full<-factor(df_counts_percat$category_full, levels=unique(df_counts_percat$category_full[order(df_counts_percat$order, decreasing=FALSE)]))
df_counts_percat$dataset <- factor(df_counts_percat$dataset, levels = c("Opentargets", "Onsides"))  # 

bp <- ggplot(df_counts_percat, aes(fill=n, y=n, x=category_full)) + 
  geom_col() +
  coord_flip() +
  facet_wrap(~dataset) +
  theme_classic() +
  scale_fill_gradient2() +
  labs(y = "Number of Distinct Drug-Side Effect Pairs",
       x = "Category" ) +theme(legend.position="none") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 3.5) 
ggsave(bp, file=paste0(plots_dir,'/Sup_Fig3_side_effect_countbycategory.png'), width=7, height=6, dpi=300)

# Supplementary Fig 4 and 5

#a. Run univariate model for each genetic feature 
lapply(c('Opentargets','Onsides'), function(dataset){
  
  datafile<-list.files(path='Data/', pattern='dataset_se.mi_withgenetics', full.names=T)
  datafile=datafile[grepl(paste0(dataset), datafile)]
  dataset1<-fread(datafile, data.table=F)
  
  univarresults_byphecodecat<-do.call(rbind,mclapply(geneticpredictors, function(pheno){
    univ_bycat<-do.call(rbind,mclapply(c(unique(dataset1$category_full)), function(category_single){
      dataset_cat=dataset1 %>% filter(category_full==category_single)
      genes_pred=dataset_cat %>% distinct(gene,.data[[pheno]]) %>% filter(.data[[pheno]]!=0)
      PT_pred=dataset_cat %>% distinct(parentterm,.data[[pheno]]) %>% filter(.data[[pheno]]!=0)
      genes_pred_outcome=dataset_cat %>% distinct(gene,.data[[pheno]],.data[[outcome]]) %>% filter(.data[[pheno]]!=0) %>% filter(.data[[outcome]]==1)
      PT_pred_outcome=dataset_cat %>% distinct(parentterm,.data[[pheno]],.data[[outcome]]) %>% filter(.data[[pheno]]!=0) %>% filter(.data[[outcome]]==1)
      PT_pred_gene_outcome=dataset_cat %>% distinct(gene,parentterm,.data[[pheno]],.data[[outcome]]) %>% filter(.data[[pheno]]!=0) %>% filter(.data[[outcome]]==1)
      
      se_pred_gene_outcome=dataset_cat %>% distinct(parentterm,.data[[pheno]],.data[[outcome]]) %>% filter(.data[[outcome]]==1)
      
      model1 =glm(as.formula(paste0(outcome, ' ~', pheno)), data=dataset_cat,family = 'binomial')
      ta1<-rbind(cbind.data.frame(outcome=outcome,Analysis=paste('Bycategory', category_single), predictors=names(model1$coefficients)[!grepl('category|Inter', names(model1$coefficients))],
                                  No.genes=length(unique(genes_pred$gene)), No.parentterms=length(unique(PT_pred$parentterm)),
                                  No.genes.outcome=length(unique(genes_pred_outcome$gene)), No.parentterms.outcome=length(unique(PT_pred_outcome$parentterm)), total.genes=length(unique(dataset_cat$gene)),  total.parentterm=length(unique(dataset_cat$parentterm)), 
                                  OR=round(exp(summary(model1)$coefficient[,1])[!grepl('category|Inter', names(model1$coefficients))],3), 
                                  lowerCI=round(exp(summary(model1)$coefficient[,1][!grepl('category|Inter', names(model1$coefficients))]-(1.96* summary(model1)$coefficient[,2][!grepl('category|Inter', names(model1$coefficients))])),3),
                                  upperCI=round(exp(summary(model1)$coefficient[,1][!grepl('category|Inter', names(model1$coefficients))]+(1.96* summary(model1)$coefficient[,2][!grepl('category|Inter', names(model1$coefficients))])),3),P.val=summary(model1)$coefficient[,4][!grepl('category|Inter', names(model1$coefficients))]))
      return(ta1)
      
    }, mc.cores=10))
    
  }, mc.cores=4))
  
  
  univar_data_outcome <- univarresults_byphecodecat  %>% mutate(upperCI=as.numeric(upperCI),Pval_sig=as.factor(ifelse(P.val>=0.05,0,1)), Analysis=gsub('Bycategory ','',Analysis))
  if(length(unique(univar_data_outcome$Pval_sig[!is.na(univar_data_outcome$Pval_sig)]))==1){
    shape_vals=c(19)
  } else{
    shape_vals=c(1,19)
  }
  
  univar_data_outcome <- FindReplace(data = univar_data_outcome, Var = "predictors", replaceData = Predictor_names, from = "predictors", to = "New", exact = TRUE)
  univar_data_outcome <-inner_join(univar_data_outcome,Category_names) %>% arrange(Category)
  univar_data_outcome$order=rep(length(unique(univar_data_outcome$Analysis)):1, each=4)
  univar_data_outcome$predictors<-factor(univar_data_outcome$predictors, levels=unique(univar_data_outcome$predictors[order(univar_data_outcome$order, decreasing=T)]))
  
  maxor=ifelse(max(univar_data_outcome$OR[!is.na(univar_data_outcome$OR)])>35, 35, max(univar_data_outcome$OR[!is.na(univar_data_outcome$OR)]))
  maxor=maxor+2
  univar_data_outcome$upperCI_cut<-ifelse(univar_data_outcome$upperCI> maxor, maxor, NA)
  
  univar_data_outcome$CI= ifelse( univar_data_outcome$upperCI>100, paste0(round(univar_data_outcome$OR,1),' (', round(univar_data_outcome$lowerCI,1), ' - ', format(round(univar_data_outcome$upperCI,1), nsmall = 1),')')
                                  , paste0(round(univar_data_outcome$OR,1),' (', round(univar_data_outcome$lowerCI,1), ' - ', round(univar_data_outcome$upperCI,1),')'))
  
  
  fp <- create_forest_plot(df = univar_data_outcome,
                           yvar = "predictors", 
                           shapevar = "Pval_sig",
                           colorvar = "Category",
                           ci_lower = "lowerCI",
                           ci_upper = "upperCI",
                           y_levels = univar_data_outcome$predictors,
                           shape_vals = shape_vals,
                           coord_limit = maxor)
  
  fp <- fp  +  facet_wrap(~Analysis) & theme(legend.position = "bottom")
  
  # if CI are too long add arrows at upper cut off 
  if (any(!is.na(univar_data_outcome$upperCI_cut))) {
    fp <- add_arrows(fp, maxor, "predictors")
  }

  if(dataset=='Opentargets'){
    ggsave(fp, file=paste0(plots_dir,'/Sup_Fig4_Forestplot_Univar_bycategory_',dataset,'.png'), width = 8, height=8, dpi=300)
  } else{
    ggsave(fp, file=paste0(plots_dir,'/Sup_Fig5_Forestplot_Univar_bycategory_',dataset,'.png'), width = 8, height=8, dpi=300)
  }
})


## Sup Fig 6 and 14, Table S1 and S3

### Sup Figure 6 and 
lapply(c('All','DOE'), function(analysis){
  dataset='Opentargets';valdataset='onside'
  Mixedmodel_all1<-do.call(rbind,lapply(c(paste0('CVsample',rep(1:5))), function(CVsample){
    if(analysis=='All'){
      mixedmodelresults=fread(paste0(scorefolder, '/Mixedeffect_weighted_model_ADR_severity_',dataset,'_outcome_',outcome_trained,'_',CVsample,'.txt'), data.table=F)
    }else {
      mixedmodelresults=fread(paste0(scorefolder, '/Mixedeffect_weighted_model_ADR_severity_',dataset,'_outcome_',outcome_trained,'_',CVsample,'_DOE_matchmechanism.txt'), data.table=F)
    }
    mixedmodelresults=mixedmodelresults[!grepl('category|Intercept|moainhibitor', mixedmodelresults$Predictor),]
    return(mixedmodelresults)
  }))
  Mixedmodel_all1$beta_sig<-as.factor(ifelse(Mixedmodel_all1$P.val>=0.05,0,1))
  Mixedmodel_all1 <- FindReplace(data = Mixedmodel_all1, Var = "Predictor", replaceData = Predictor_names,
                               from = "predictors", to = "New", exact = TRUE)

  Mixedmodel_all1<-inner_join(Mixedmodel_all1, Category_names, by=c('Predictor'='predictors'))
  Mixedmodel_all<- Mixedmodel_all1%>% arrange(Category, Predictor,CV)
  Mixedmodel_all$Predictor<- factor(Mixedmodel_all$Predictor, levels=unique(Mixedmodel_all$Predictor))
  Mixedmodel_all$CV<-as.factor(Mixedmodel_all$CV)
  
  if(length(unique(Mixedmodel_all$beta_sig[!is.na(Mixedmodel_all$beta_sig)]))==1){
    shape_vals=c(19)
  } else{
    shape_vals=c(1,19)
  }
  minstart=ifelse(min(Mixedmodel_all$lowerCI)>0,-0.1,min(Mixedmodel_all$lowerCI)>0)
  # Create forest plot 
  fp <- Mixedmodel_all %>% as_tibble() %>%
    ggplot(aes(y=Predictor, x=beta,color=CV))+
    geom_point(size=1.4, stroke=0.5,aes(color=CV, shape =beta_sig), position =position_dodge(0.5))  +
    scale_shape_manual(values=shape_vals) +guides(shape = "none") + 
    geom_linerange(aes(xmin=lowerCI, xmax=upperCI,color=CV),position = position_dodge(0.5)) +
    scale_y_discrete(limits = rev(levels(Mixedmodel_all$Predictor))) +xlim(minstart, max(Mixedmodel_all$upperCI)) +
    xlab('Beta coefficients (95% CI)') + ylab('') +#  ggtitle(paste(dataset, outcome,'\n', title, drugtype,'\n', take )) +
    geom_vline(xintercept=0, linetype='longdash', color='red') +
    theme(axis.text=element_text(size=13,color = "black"), axis.title=element_text(size=14,color = "black"),
          axis.line = element_line(colour = 'black', size = 0.5), 
          axis.ticks = element_line(colour = "black", size = 0.5),
          axis.ticks.length = unit(0.25, "cm"),
          plot.margin=margin(10,20,10,10)) +
    theme_bw()
  
  if(analysis=='All'){
    filename=paste0(plots_dir,'/Sup_Fig6_Mixedmodel_',dataset,'_outcometrained_',outcome_trained,'_multivariatemodel_results.png')  
  } else{
    filename=paste0(plots_dir,'/Sup_Fig13_Mixedmodel_',dataset,'_outcometrained_',outcome_trained,'_multivariatemodel_results_DOE.png')  
  }
  ggsave(fp, file= filename, width = 7, height=6, dpi=300)

})


## Sup Figure 7


#### Violin and Bar plots
geneticpredictors2=c( 'clinicalvariant_1','clinicalvariant_2','clinicalvariant_3','gwastrait','geneburden','singlevar')
Predictor_names2= data.frame(from = geneticpredictors2, to = c("Clinical Variant 1 predictor","Clinical Variant 2 predictors","Clinical Variant 3 predictors","GWA trait", "Gene Burden","Single Variant"))


lapply(c('Onsides','Opentargets'), function(dataset){
  
  lapply(c('All','DOE_inhibitor','DOE_activator'), function(analysis){
    print(paste(analysis ))
    
    if(analysis=='DOE_inhibitor'){
      
      
      Genescores_beta_genetic_CV=fread(paste0(scorefolder, 'All_genescoresum_across_all_predictors_',dataset,'_outcome_',outcome_trained,'_sideeffectproject_DOE.txt.gz'), data.table=F)
      Genescores_beta_genetic_CV$genescoresum_pos= ifelse(Genescores_beta_genetic_CV$genescoresum >0,  Genescores_beta_genetic_CV$genescoresum , 0)
      Dataset_genescores2= Genescores_beta_genetic_CV%>% select(-genescoresum) %>% dplyr::rename(genescoresum=genescoresum_pos) %>% arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% 
        mutate(percent=order/length(genescoresum) *100) %>% dplyr::rename(gwastrait=gwastrait_doe,clinicalvariant = clinicalvariant_DOE,geneburden = geneburden_doe, singlevar = singlevar_doe)
      Dataset_genescores1= Dataset_genescores2%>% 
        dplyr::mutate(across(all_of(geneticpredictors), ~ ifelse(. < 0, 0, .))) # restrict to just LOF/positive predictors 
      
    }else if(analysis=='DOE_activator'){
      
      Genescores_beta_genetic_CV=fread(paste0(scorefolder, 'All_genescoresum_across_all_predictors_',dataset,'_outcome_',outcome_trained,'_sideeffectproject_DOE.txt.gz'), data.table=F)
      Genescores_beta_genetic_CV$genescoresum_neg= ifelse(Genescores_beta_genetic_CV$genescoresum < 0,  Genescores_beta_genetic_CV$genescoresum , 0)
      Dataset_genescores2= Genescores_beta_genetic_CV%>% select(-genescoresum) %>% dplyr::rename(genescoresum=genescoresum_neg) %>% mutate(genescoresum=genescoresum*-1) %>% arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>%
        mutate(percent=order/length(genescoresum) *100) %>%   dplyr::rename(gwastrait=gwastrait_doe, clinicalvariant = clinicalvariant_DOE,geneburden = geneburden_doe, singlevar = singlevar_doe)
      Dataset_genescores1= Dataset_genescores2%>% 
        dplyr::mutate(across(all_of(geneticpredictors), ~ ifelse(. > 0, 0, -.)))# restrict to GOF/negative predictors 
    } else { # For All analysis
      Genescores_beta_genetic_CV=fread(paste0(scorefolder, 'All_genescoresum_across_all_predictors_',dataset,'_outcome_',outcome_trained,'_sideeffectproject.txt.gz'), data.table=F)
      Dataset_genescores1= Genescores_beta_genetic_CV %>%arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% mutate(percent=order/length(genescoresum) *100)
    }
    #calculate number of clinical variant pred that contributed to CV predictor
    Dataset_genescores1$no.cv=0
    CV_beta_buffer= min(Dataset_genescores1$clinicalvariant[Dataset_genescores1$clinicalvariant!=0]) +0.15 #cv beta plus .15 for varaibility between cv samples. X2 to get 2 predictors etc.
    Dataset_genescores1$no.cv<-ifelse(Dataset_genescores1$clinicalvariant>0 & Dataset_genescores1$clinicalvariant <CV_beta_buffer, 1,   Dataset_genescores1$no.cv) #
    Dataset_genescores1$no.cv<-ifelse(Dataset_genescores1$clinicalvariant>CV_beta_buffer & Dataset_genescores1$clinicalvariant <CV_beta_buffer *2, 2,   Dataset_genescores1$no.cv)
    Dataset_genescores1$no.cv<-ifelse( Dataset_genescores1$clinicalvariant >CV_beta_buffer *2, 3,   Dataset_genescores1$no.cv)
    Dataset_genescores1$no.cv_additional=ifelse( Dataset_genescores1$no.cv!=0, Dataset_genescores1$no.cv-1,  Dataset_genescores1$no.cv)
    ## calculate no predictors  
    Dataset_genescores1$no.pred=rowSums(Dataset_genescores1[,geneticpredictors] != 0)  
    Dataset_genescores1$no.pred=Dataset_genescores1$no.pred+ Dataset_genescores1$no.cv_additional
    Dataset_genescores1$clinicalvariant_1<-ifelse(Dataset_genescores1$no.cv==1,Dataset_genescores1$clinicalvariant,0)
    Dataset_genescores1$clinicalvariant_2<-ifelse(Dataset_genescores1$no.cv==2,Dataset_genescores1$clinicalvariant,0)
    Dataset_genescores1$clinicalvariant_3<-ifelse(Dataset_genescores1$no.cv==3,Dataset_genescores1$clinicalvariant,0)
    Dataset_genescores1_df=Dataset_genescores1 %>% select(contains(c(geneticpredictors))) %>% select(-clinicalvariant)
    Dataset_genescores1_df1=cbind(percent=Dataset_genescores1$percent, Dataset_genescores1_df) 
    Dataset_genescores1_df2=pivot_longer(data = Dataset_genescores1_df1, cols = -c('percent'), names_to = "Predictor", values_to = "weights")
    Dataset_genescores1_df2_1=Dataset_genescores1_df2[!is.na(Dataset_genescores1_df2$weights),]
    Dataset_genescores1_df3=Dataset_genescores1_df2_1[Dataset_genescores1_df2_1$weights!=0,]  
    #2. Make plot
    Dataset_genescores1_df3=as.data.frame(Dataset_genescores1_df3)
    Allfiletypes1 <- FindReplace(data = Dataset_genescores1_df3, Var = "Predictor", replaceData = Predictor_names2,
                                 from = "from", to = "to", exact = TRUE)
    Allfiletypes1$Predictor<-factor(Allfiletypes1$Predictor, levels=unique(Allfiletypes1$Predictor[order(Allfiletypes1$weights, decreasing=TRUE)]))
    weights_an= Allfiletypes1 %>% distinct(Predictor, weights) %>% group_by(Predictor) %>% dplyr::summarise(median.weights=median(weights))%>% arrange(desc(median.weights)) 
    print(weights_an)
    Weights_all_order=weights_an%>% arrange(desc(median.weights)) 
    Weights_all_order$order=(nrow(Weights_all_order):1)+0.2
    Allfiletypes1$Predictor<- factor(Allfiletypes1$Predictor, levels = unique(Weights_all_order$Predictor[order(Weights_all_order$median.weights, decreasing=T)]))
    ## add sample size to plot 
    sample_size_allpred = Allfiletypes1%>% group_by(Predictor) %>% tally() %>% dplyr::rename(num=n) 
    Allfiletypes12=Allfiletypes1 %>% left_join(sample_size_allpred) %>%
      mutate(myaxis = paste0(Predictor, "\n", "(n=", num,")"))
    Weights_all_order1=Weights_all_order %>% left_join(sample_size_allpred) %>%
      mutate(myaxis = paste0(Predictor, "\n", "(n=", num,")"))
    Allfiletypes12$myaxis<- factor(Allfiletypes12$myaxis, levels = unique(Weights_all_order1$myaxis[order(Weights_all_order1$median.weights, decreasing=T)]))
    
    if(analysis=='All') {
      xlabel=('Percentile of SE-GPS (%)')
    } else{
      xlabel=('Percentile of SE-GPS-DOE (%)')
    }
    
    plot3<- ggplot(data = Allfiletypes12, mapping = aes(x =percent,y = myaxis, fill=Predictor)) +
      geom_violin(position = position_dodge(width = 0.2), scale = "width", width=0.3,adjust = 0.7, trim = T, size=0.2)+
      scale_x_continuous( limits=c( min(Allfiletypes12$percent), 100)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "point", shape=21, size=3.5, color='black') +
      theme(legend.position="none") +   guides(fill = F, color = F) + scale_y_discrete(limits = rev(levels(Allfiletypes12$myaxis))) +
      geom_text(data= Weights_all_order, aes(x=min(Allfiletypes12$percent), y=order) , label = paste(round(Weights_all_order1$median.weights,3)),size=3) +
      xlab(xlabel) + ylab('')  + 
      theme(axis.text=element_text(size=14), axis.title=element_text(size=15) )  +
      scale_fill_manual(values=c('Clinical Variant 1 predictor'= '#DB72FB','Clinical Variant 2 predictors'= '#00BA38','Clinical Variant 3 predictors'= '#00B9E3', 'Single Variant'='#619CFF', 'Gene Burden'= '#00C19F',
                                 'GWA trait'='#F8766D'))+  theme_classic() +
      theme(axis.line = element_line(colour = 'black', size = 0.5), 
            axis.ticks.y = element_line(colour = "black", size = 0.5),
            axis.ticks.x = element_line(colour = "black", size = 0.5),
            axis.ticks.length = unit(0.25, "cm")) +
      theme(axis.text.x=element_text(color = "black",),
            axis.text.y=element_text(color = "black"),
            axis.title=element_text(color = "black")) +
      theme(plot.margin=margin(10,20,10,10)) 
    
    
    if(analysis=='All' & dataset=='Opentargets') {
      ggsave(plot3 , file=paste0(plots_dir, '/Sup_Fig7_Violin_plot_',dataset,'_',analysis,'.png'), width = 6, height=6, dpi=300)
    } else if(analysis=='All' & dataset=='Onsides'){
      ggsave(plot3 , file=paste0(plots_dir, '/Sup_Fig9_Violin_plot_',dataset,'_',analysis,'.png'), width = 6, height=6, dpi=300)
    }
    
    
    ### bar plots
    
    print(head(Dataset_genescores1))
    
    if(analysis=='DOE'){
      binmax=1.8
    } else{
      binmax=2.1}
    binsize =0.3
    breaks1=seq(0,binmax,binsize)
    labels1=paste(breaks1,'-', breaks1+binsize)
    breaks1=c(breaks1, Inf)
    labels1=gsub(paste0(binmax,' - ', binmax+binsize),paste0(binmax,'+'), labels1)
    Dataset_genescores1_group=Dataset_genescores1 %>% filter(genescoresum!=0)
    Dataset_genescores1_group$group <- cut(Dataset_genescores1_group$genescoresum, breaks = breaks1, labels = labels1, right = FALSE, include.lowest=F)
    
    Dataset_genescores1_df=Dataset_genescores1_group %>% select(contains(c(geneticpredictors))) %>% select(-clinicalvariant)
    Dataset_genescores1_df1=cbind(group=Dataset_genescores1_group$group, no.pred=Dataset_genescores1_group$no.pred, percent= Dataset_genescores1_group$percent, Dataset_genescores1_df)
    
    Dataset_genescores1_df2=data.frame(pivot_longer(data = Dataset_genescores1_df1, 
                                                    cols = -c(1:3), 
                                                    names_to = "Predictor", 
                                                    values_to = "weights"))
    
    Dataset_genescores1_df2_count1=data.frame(Dataset_genescores1_df2 %>% group_by(group,no.pred,percent) %>% filter(weights!=0) %>%  dplyr::summarise(Predictor=paste0(Predictor, collapse=','))) 
    Dataset_genescores1_df2_count2=data.frame(Dataset_genescores1_df2_count1 %>% group_by(group,no.pred,Predictor) %>%  tally() %>%   dplyr::mutate(numbering = row_number()) %>%
                                                separate_rows(Predictor,sep=',') %>% dplyr::rename(value=n))
    
    Dataset_genescores1_df2_count2 <- FindReplace(data = Dataset_genescores1_df2_count2, Var = "Predictor", replaceData = Predictor_names2,
                                                  from = "from", to = "to", exact = TRUE)
    Pred_bygroup <- tibble(Predictor = rep(unique(Dataset_genescores1_df2_count2$Predictor),length(unique(as.character(Dataset_genescores1_df2_count2$group)))),
                           group =  rep(unique(as.character(Dataset_genescores1_df2_count2$group)),length(unique(as.character(Dataset_genescores1_df2_count2$Predictor)))))
    Dataset_genescores1_df2_count2_pred<-full_join(Pred_bygroup, Dataset_genescores1_df2_count2) %>% arrange(group,Predictor )
    Dataset_genescores1_df2_count2$group<-as.character(Dataset_genescores1_df2_count2$group)
    
    
    if(analysis=='All') {
      xlabel=('Number of genetic features contributing to the SE-GPS')
    } else{
      xlabel=('Number of genetic features contributing to the SE-GPS-DOE')
    }
    
    bp <- ggplot(Dataset_genescores1_df2_count2, aes(fill=Predictor, y=value, x=no.pred)) + 
      geom_bar(position="stack", stat="identity") + 
      geom_hline(yintercept=0, color = "grey") +
      facet_wrap(~group,scales="free_y") +
      theme_classic() + xlab(xlabel) + ylab('Count') +labs(fill = "Feature") +    
      scale_fill_manual(values=c('Clinical Variant 1 predictor'= '#DB72FB','Clinical Variant 2 predictors'= '#00BA38','Clinical Variant 3 predictors'= '#00B9E3', 'Single Variant'='#619CFF', 'Gene Burden'= '#00C19F',
                                 'GWA trait'='#F8766D'))
    
    if(analysis=='All' & dataset=='Opentargets') {
      ggsave(bp, file=paste0(plots_dir,'/Sup_Fig7_Barplot_counts_genescoresumbins_',dataset,'_',analysis,'.png'), width=7, height=6, dpi=300)
    }else if(analysis=='All' & dataset=='Onsides') {
      ggsave(bp, file=paste0(plots_dir,'/Sup_Fig9_Barplot_counts_genescoresumbins_',dataset,'_',analysis,'.png'), width=7, height=6, dpi=300)
    }else if(analysis != 'All' & dataset=='Opentargets'){
      ggsave(bp, file=paste0(plots_dir,'/Sup_Fig14_barplot_counts_genescoresumbins_',dataset,'_',analysis,'.png'), width=7, height=6, dpi=300)
    }else if(analysis != 'All' & dataset=='Onsides'){
      ggsave(bp, file=paste0(plots_dir,'/Sup_Fig15_barplot_counts_genescoresumbins_',dataset,'_',analysis,'.png'), width=7, height=6, dpi=300)
    }
    
  })
})

### Supplementary Fig. 8 
dataset='Opentargets'
OTGenescores_beta_genetic_CV=fread(paste0(scorefolder, 'All_genescoresum_across_all_predictors_',dataset,'_outcome_senomi_phase4_sideeffectproject.txt.gz'), data.table=F)
OTGenescores_beta_genetic_CV_genept=OTGenescores_beta_genetic_CV %>% distinct(gene, parentterm, genescoresum) 
OTGenescores_beta_genetic_CV_geneptnonzero=OTGenescores_beta_genetic_CV %>% distinct(gene, parentterm, genescoresum) %>% filter(genescoresum>0)

dist=ggplot(OTGenescores_beta_genetic_CV_genept, aes(x = genescoresum)) +
  geom_density(color = "blue", fill = "lightblue", alpha = 0.5) +
  theme_classic() +
  labs(title = "Distribution of SE-GPS", x = "SE-GPS", y = "Density")
ggsave(dist, file=paste0(plots_dir,'/Sup_Fig8a_distribution_scores_density_',dataset,'.png'), width = 6.5, height=5, dpi=300)

dist=ggplot(OTGenescores_beta_genetic_CV_geneptnonzero, aes(x = genescoresum)) +
  geom_density(color = "blue", fill = "lightblue", alpha = 0.5) +
  theme_classic() +
  labs(title = "Distribution of SE-GPS", x = "SE-GPS", y = "Density")
ggsave(dist, file=paste0(plots_dir,'/Sup_Fig8b_distribution_scores_density_nonzero',dataset,'.png'), width = 6.5, height=5, dpi=300)


### Supplementary Fig 18

library(VennDiagram)

ot=fread('Data/Opentargets_dataset_se.mi_withgeneticsfiltered_5_All_parentterm_drugse_final.txt.gz', data.table=F)
onsides=fread('Data/Onsides_dataset_se.mi_withgeneticsfiltered_5_All_parentterm_drugse_final.txt.gz', data.table=F)
ot$gene_parentterm=paste(ot$gene, '_',ot$parentterm)
onsides$gene_parentterm=paste(onsides$gene, '_',onsides$parentterm)

# overlap data including side effects

ot_se=ot %>% filter(senomi_phase4==1) 
onsides_se=onsides %>% filter(senomi_phase4==1) 

overlap_data_se = list(`Open Targets` = unique(ot_se$gene_parentterm), `OnSIDES`=unique(onsides_se$gene_parentterm))

venn.diagram(
  x = overlap_data_se,
  category.names = c("Open Targets", "OnSIDES"),
  fill = c("blue", "green"),  
  alpha = 0.30,
  print.mode = c("raw", "percent"),
  cat.cex = 1,
  cex = 1,
  fontfamily = "sans",
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily='sans',
  cat.dist = c(0.05, 0.05),  # <-- Push labels out more
  margin = 0.07,  # <-- Increase to give more space
  filename = "results/plots/Sup_Fig18_venn_gene-parenttermoverap_sideeffects.png",
  imagetype = "png",
  output = TRUE
)


end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

