
# Load required libraries
library(data.table)
library(parallel)
library(dplyr)
library(tidyr)

# Set required paths

scorefolder<-'results/construct_scores/'
tables_dir='results/tables/'
dataset='Opentargets'; outcome='senomi_phase4'


if (!dir.exists(scorefolder)) {
  dir.create(scorefolder, recursive = TRUE)
}

if (!dir.exists(tables_dir)) {
  dir.create(tables_dir, recursive = TRUE)
}

start.time <- Sys.time()
print(start.time)

set.seed(18)

geneticpredictors= c('clinicalvariant', 'gwastrait','geneburden','singlevar') 

#######################################
# 1. Create 5 fold CV training sets 
#######################################

 # Open Target drug-genetic dataset
datafile<-list.files(path='Data/', pattern='Opentargets_dataset_se.mi_withgenetics', full.names=T)
opentarget_fulldataset<-fread(datafile, data.table=F)
opentarget_fulldataset$gene_phenotype=paste0(opentarget_fulldataset$gene,'-',opentarget_fulldataset$parentterm)

# Create five equally sized groups of randomly assigned unique gene-phenotype pairs -  5 cross validation 20% sets
# Each group is test set and training set will be remaining 80%
uniquegene_phenotype1=opentarget_fulldataset%>% distinct(gene_phenotype)
uniquegene_phenotype_ss1 <- sample(rep(1:5, diff(floor(nrow(uniquegene_phenotype1) * c(0, 0.2, 0.4, 0.6, 0.8, 1)))))
splitdata_gp1 <- split(uniquegene_phenotype1,uniquegene_phenotype_ss1)

# create cross validation test set
for(i in 1:length(splitdata_gp1)){
  data_split=splitdata_gp1[[i]]
  test_set <- opentarget_fulldataset[opentarget_fulldataset$gene_phenotype  %in% data_split$gene_phenotype, ]
  filename=paste0('results/construct_scores/',dataset,'_se_dataset_sampled20_CVsample',i,'.txt.gz')
  print(filename)
  write.table(test_set, gzfile(filename), sep='\t', quote=F, row.names=F)
}


# unload and reload required libraries
unloadNamespace("tidyr")
unloadNamespace("data.table")
unloadNamespace("parallel")
unloadNamespace("dplyr")
unloadNamespace("tidyr")

library("Matrix")
library("Rcpp")
library("StanHeaders")
library("rstan")
options(mc.cores = parallel::detectCores())
library("lme4")
library("brms")
library(data.table)


#######################################
# 2. Run mixed model in 80% training set 
#######################################

opentarget_fulldataset<-fread(datafile, data.table=F)
opentarget_fulldataset$ID=paste0(opentarget_fulldataset$drugname,'-',opentarget_fulldataset$gene,'-',opentarget_fulldataset$parentterm)

# Combine the mixed effects results model across 5 cross-validation (CV) samples
Mixedmodel_all<-do.call(rbind, lapply(c(paste0('CVsample',rep(1:5))), function(CVsample){
  CVdataset=fread(paste0(scorefolder,dataset,'_se_dataset_sampled20_',CVsample,'.txt.gz'),data.table=F) #Load the 20% test Open target dataset for cross validated sample
  CVdataset$ID=paste0(CVdataset$drugname,'-',CVdataset$gene,'-',CVdataset$parentterm)   # Create a unique identifier for each drug-gene-phenotype triplet
  dataset_sampled80=subset(opentarget_fulldataset, !(ID %in% CVdataset$ID)) # Generate the remaining 80% training dataset by excluding the current test set entries
  # Load the ADR severity scores
  ADR_severity1=fread(paste0('Data/Adr_severity_scores_phecodeX.txt'), data.table=F) 
  dataset_sampled80_ADR=merge(dataset_sampled80, ADR_severity1, all.x=T)
  dataset_sampled80_ADR$genepheno=paste0(dataset_sampled80_ADR$gene,'-', dataset_sampled80_ADR$parentterm)
  # Round severity scores (scaled by 100) for weighting the regression
  dataset_sampled80_ADR$severity_score_round=round(dataset_sampled80_ADR$severity_score*100)
  Predictors=paste0(paste(geneticpredictors,collapse='+'), '+category ')
  #Run mixedmodel regression across each CV training dataset to get weights with drugname as a random effect
  mod_se_mixed <- glmer(as.formula(paste0(outcome, " ~ ", paste(Predictors), " +(1 | drugname)")),
                        family = "binomial",
                        data = dataset_sampled80_ADR,
                        control = glmerControl(optimizer = "bobyqa"),nAGQ=0L,
                        weights=severity_score_round)
  
  # Extract regression results (coefficients, confidence intervals, p-values)
  mod_se_mixed_results1 <- cbind.data.frame(beta = summary(mod_se_mixed)$coefficient[, c(1)],
                                            lowerCI = summary(mod_se_mixed)$coefficient[, c(1)] - 1.96 * summary(mod_se_mixed)$coefficient[, c(2)],
                                            upperCI = summary(mod_se_mixed)$coefficient[, c(1)] + 1.96 * summary(mod_se_mixed)$coefficient[, c(2)],
                                            P.val = summary(mod_se_mixed)$coefficient[, c(4)])
  
  mod_se_mixed_results1$Predictor=rownames(mod_se_mixed_results1)
  mod_se_mixed_results1$CV=CVsample
  
  write.table(mod_se_mixed_results1, paste0(scorefolder,'/Mixedeffect_weighted_model_ADR_severity_',dataset,'_outcome_',outcome,'_',CVsample,'.txt'), sep='\t',quote=F)
  return(mod_se_mixed_results1)
}))

print(paste('load all packges...'))
print(paste(head(Mixedmodel_all)))

# load required libraries
library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(DataCombine)
library(logistf)
library(parallel)


Sup_table_weights<-Mixedmodel_all %>% dplyr::mutate(CI=paste0(round(lowerCI,2), ' - ', round(as.numeric(upperCI),2)), beta=round(beta,2)) %>%
  select(CV, Predictor, beta, CI, P.val)  %>% dplyr::rename(`P-value`=`P.val`) %>% arrange(CV)

write.table(Sup_table_weights, paste0(tables_dir,'/Table_S1_Mixedmodel_',dataset,'_multivariatemodel_results.txt') , sep='\t', quote=F, row.names=F)



##############################################################################
# 3. Create scores - apply betas from mixed model as weights to CV test set
##############################################################################

opentarget_fulldataset<-fread(datafile, data.table=F)
opentarget_fulldataset$gene_phenotype=paste0(opentarget_fulldataset$gene,'-',opentarget_fulldataset$parentterm)

lapply(c(paste0('CVsample',rep(1:5))), function(CVsample){
  
  mixedmodelresults=fread(paste0(scorefolder, '/Mixedeffect_weighted_model_ADR_severity_',dataset,'_outcome_',outcome,'_',CVsample,'.txt'), data.table=F)
  mixedmodelresults=mixedmodelresults[!grepl('category|Intercept', mixedmodelresults$Predictor),] # extract coefficient for each genetic predictor as weights
  #Use weights from Mixedmodel and create score using remaining 20% of data for each CV
  Model_weights=as.data.frame(t(as.numeric(mixedmodelresults$beta)))
  names(Model_weights)=mixedmodelresults$Predictor
  #extract gene, parentterm and genetic predictor columns
  CVdataset=fread(paste0(scorefolder,dataset,'_se_dataset_sampled20_',CVsample,'.txt.gz'),data.table=F) #20% test Open target dataset
  CVdataset_gene_phenotype=CVdataset[, grepl(paste0('\\bgene\\b|parentterm|',paste0(geneticpredictors, collapse='\\b|')), colnames(CVdataset))] %>% distinct() 
  # Multiply each predictor value in the test set by its corresponding beta coefficient weight
  Genescores_beta_genetic_values=mapply("*", CVdataset_gene_phenotype[intersect(names(CVdataset_gene_phenotype), names(Model_weights))],
                                        Model_weights[intersect(names(CVdataset_gene_phenotype), names(Model_weights))])
  gene_parentterm=CVdataset_gene_phenotype %>% select(gene, parentterm)
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #Construct SE-GPS as the sum across all predictor values for each gene-phenotype pair
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm")])
  write.table(Genescores_beta_genetic, paste0(scorefolder, 'Genescore_sum_across_predictor_',dataset,'_outcome_',outcome,'_',CVsample,'.txt'), sep='\t', row.names=F, quote=F )
})

## Choose the CV sample with the max OR for validation set
max_filetype<-do.call(rbind,lapply(c(paste0('CVsample',seq(1:5))), function(CVsample){   
  CVdataset=fread(paste0(scorefolder,dataset,'_se_dataset_sampled20_',CVsample,'.txt.gz'),data.table=F) #20% test Open target dataset
  genescorefile=fread(paste0(scorefolder, 'Genescore_sum_across_predictor_',dataset,'_outcome_',outcome,'_',CVsample,'.txt'), data.table=F)
  # combine genescore sum file with drug data for each test set 
  Dataset_genescores=inner_join(CVdataset[c('drugname','gene','parentterm','category','category_full','mi','se','senomi_phase4')] ,genescorefile, by=c('gene','parentterm') )
  model1 =glm(as.formula(paste0(outcome, ' ~ genescoresum + category')), data=Dataset_genescores,family = 'binomial')
  mod_output<-rbind(cbind.data.frame(CV=CVsample,outcome=outcome,OR=exp(summary(model1)$coefficient[2,1]),lowerCI=exp(summary(model1)$coefficient[2,1]-(1.96* summary(model1)$coefficient[2,2])),upperCI=exp(summary(model1)$coefficient[2,1]+(1.96* summary(model1)$coefficient[2,2])),P.val=summary(model1)$coefficient[2,4]))
  return(mod_output)
}))

## take CV with the max OR
Sup_table<-max_filetype %>% mutate(CI=paste0(round(lowerCI,2), ' - ', round(as.numeric(upperCI),2)), OR=round(OR,2)) %>%
  select(CV, OR, CI, P.val)  %>% dplyr::rename(`P-value`=`P.val`) %>% arrange(desc(OR))

write.table(Sup_table, paste0(tables_dir, 'Table_S2_associations_across_CVtest_sets.txt'), sep='\t', quote=F, row.names=F)

## For each CV-fold, merge the SE-GPS with the corresponding CV 20% drug dataset at the gene-phenotype level to combine the SE-GPS with drug side effect information
## combine into a single data frame (100% of data)
Combine_genescores<-lapply(paste0('CVsample', seq(1:5)), function(CVsample){
  CVdataset=fread(paste0(scorefolder,dataset,'_se_dataset_sampled20_',CVsample,'.txt.gz'),data.table=F) #20% test Open target dataset
  genescorefile=fread(paste0(scorefolder, 'Genescore_sum_across_predictor_',dataset,'_outcome_',outcome,'_',CVsample,'.txt'), data.table=F)
  Genescore_sumfile_20_drugs<-inner_join(CVdataset[c('drugname','gene','parentterm','category','category_full','mi','se','senomi_phase4')], genescorefile, by = c('gene' = 'gene', 'parentterm' = 'parentterm')) %>% distinct()
  Genescore_sumfile_20_drugs$CV=CVsample
  return(Genescore_sumfile_20_drugs)
})

Combine_genescores1<-do.call(rbind.fill,Combine_genescores)
write.table(Combine_genescores1, gzfile(paste0(scorefolder, 'All_genescoresum_across_all_predictors_',dataset,'_outcome_',outcome,'_sideeffectproject.txt.gz')), sep='\t',quote=F,row.names=F)

##############################################################################
# 4. Apply method to OnSIDES dataset 
##############################################################################

lapply(c('Onsides'), function(valdataset){
  #Calculate SE-GPS in OnSIDES dataset using the OT model weights that produced the highest OR 
  OT_CV=fread(paste0(tables_dir, 'Table_S2_associations_across_CVtest_sets.txt'), data.table=F )
  samplecv=head(OT_CV[1],1)[[1]]
  # Load mixed model coefficients (predictor weights) from Open Targets for the selected CV fold
  mixedmodelresults=fread(paste0(scorefolder, '/Mixedeffect_weighted_model_ADR_severity_Opentargets_outcome_',outcome,'_',samplecv,'.txt'), data.table=F)
  mixedmodelresults=mixedmodelresults[!grepl('category|Intercept', mixedmodelresults$Predictor),]
  mixedmodel=mixedmodelresults %>% select(Predictor, beta) 
  # Load the OnSIDES validation dataset
  valdataset_dir<-list.files(path=paste0('Data/'), pattern='Onsides_dataset_se.mi_withgenetics', full.names=T)
  Validation_dataset=fread(paste0(valdataset_dir),data.table=F) 
  # Extract beta coefficients 
  Model_weights=as.data.frame(t(mixedmodel$beta))
  names(Model_weights)=mixedmodel$Predictor
  
  # Select unique gene-phenotype predictor combinations from the validation dataset
  Validation_dataset_gene_phenotype=Validation_dataset[, grepl(paste0('\\bgene\\b|parentterm|',paste0(geneticpredictors, collapse='\\b|')), colnames(Validation_dataset))] %>% distinct() 
  # Multiply each predictor value in the validation set by its corresponding OT beta coefficient weight
  gene_parentterm=Validation_dataset %>% distinct(gene, parentterm)
  Genescores_beta_genetic_values=mapply("*", Validation_dataset_gene_phenotype[intersect(names(Validation_dataset_gene_phenotype), names(Model_weights))],
                                        Model_weights[intersect(names(Validation_dataset_gene_phenotype), names(Model_weights))])
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #Construct SE-GPS as the sum across all predictor values for each gene-phenotype pair
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm")])
  # Join SE-GPS with drug-data and write full output to file
  Genescores_beta_genetic_drugs=inner_join(Validation_dataset[c('drugname','gene','parentterm','category','category_full','mi','se','senomi_phase4')],Genescores_beta_genetic, by=c('gene', 'parentterm'))
  write.table(Genescores_beta_genetic_drugs, paste0(scorefolder, 'All_genescoresum_across_all_predictors_',valdataset,'_outcome_',outcome,'_sideeffectproject.txt.gz'), sep='\t', row.names=F, quote=F )
  
  
})
print(paste('ONSIDE done'))

##############################################################################
# 5. Apply method to All genes  
##############################################################################

lapply(c('Allgenes'), function(valdataset){
  #Calculate SE-GPS in all genes dataset (19400 + genes) using the OT weights which gave the max OR 
  print(valdataset)
  OT_CV=fread(paste0(tables_dir, 'Table_S2_associations_across_CVtest_sets.txt'), data.table=F )
  samplecv=head(OT_CV[1],1)[[1]]
  # Load mixed model coefficients (predictor weights) from Open Targets for the selected CV fold
  mixedmodelresults=fread(paste0(scorefolder, '/Mixedeffect_weighted_model_ADR_severity_Opentargets_outcome_',outcome,'_',samplecv,'.txt'), data.table=F)
  mixedmodelresults=mixedmodelresults[!grepl('category|Intercept', mixedmodelresults$Predictor),]
  mixedmodelresults=mixedmodelresults[!grepl('category|Intercept',mixedmodelresults$V1),]
  # Extract beta coefficients 
  mixedmodel=mixedmodelresults %>% select(Predictor, beta) 
  print(mixedmodel)
  #All genetic features across all genes
  phenotype_allgene_pt2<-fread(paste0('Data/Allgenes_dataset_parentterm_predictors_collapsed_sideeffect_project.txt.gz'), data.table=F)
  phecodef=read.csv('Data/phecodeX_info.csv')
  phecodef = phecodef %>% dplyr::rename(category_full=category) %>% mutate(category=gsub('_.*', '',phecode), phecode=gsub('.*_', '',phecode)) %>%
    mutate(parentterm=as.character(trunc(as.numeric(phecode)))) %>%  select(parentterm, category_full, category ) %>% distinct() 
  phenotype_allgene_pt2= phenotype_allgene_pt2 %>% filter(!category_full %in% c('Neonatal','Pregnancy'))
  
  #Create combined features
  print(paste('Create combined features'))
  clincicalpred=c('clinvar','omim','hgmd')
  phenotype_allgene_pt2$clinicalvariant=0
  phenotype_allgene_pt2$clinicalvariant=0
  phenotype_allgene_pt2$clinicalvariant[rowSums(phenotype_allgene_pt2[, clincicalpred]) == 1]<-1
  phenotype_allgene_pt2$clinicalvariant[rowSums(phenotype_allgene_pt2[, clincicalpred]) == 2]<-2
  phenotype_allgene_pt2$clinicalvariant[rowSums(phenotype_allgene_pt2[, clincicalpred]) == 3]<-3
  geneburden=c('OTgb','ravargb')
  phenotype_allgene_pt2$geneburden=0
  phenotype_allgene_pt2$geneburden[rowSums(phenotype_allgene_pt2[, geneburden]) != 0]<-1
  singlevar=c('ravarsv','genebasssv')
  phenotype_allgene_pt2$singlevar=0
  phenotype_allgene_pt2$singlevar[rowSums(phenotype_allgene_pt2[, singlevar]) != 0]<-1
  gwastrait=c('l2g','eqtlpheno')
  phenotype_allgene_pt2$gwastrait=0
  phenotype_allgene_pt2$gwastrait[rowSums(phenotype_allgene_pt2[, gwastrait]) != 0]<-1
  
  Model_weights=as.data.frame(t(mixedmodel$beta))
  names(Model_weights)=mixedmodel$Predictor
  Validation_dataset_gene_phenotype=phenotype_allgene_pt2[, grepl(paste0('\\bgene\\b|parentterm|',paste0(geneticpredictors, collapse='\\b|')), colnames(phenotype_allgene_pt2))] %>% distinct() 
  # Multiply each predictor value in the All genes data set by its corresponding OT beta coefficient weight
  gene_parentterm=Validation_dataset_gene_phenotype %>% distinct(gene, parentterm)
  Genescores_beta_genetic_values=mapply("*", Validation_dataset_gene_phenotype[intersect(names(Validation_dataset_gene_phenotype), names(Model_weights))],
                                        Model_weights[intersect(names(Validation_dataset_gene_phenotype), names(Model_weights))])
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #Construct SE-GPS as the sum across all predictor values for each gene-phenotype pair
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm")])
  Genescores_beta_genetic$parentterm=as.character(Genescores_beta_genetic$parentterm)
  Genescores_beta_genetic_drugs=inner_join(phecodef,Genescores_beta_genetic, by=c( 'parentterm'))
  write.table(Genescores_beta_genetic_drugs, paste0(scorefolder, 'All_genescoresum_across_all_predictors_',valdataset,'_outcome_',outcome,'_sideeffectproject.txt.gz'), sep='\t', row.names=F, quote=F )
  print(paste(valdataset, 'done'))
  
})

print('withdrawn dataset')

##############################################################################
# 6. Create a severe drug dataset using boxed warnings and withdrawn drugs for
#Open Targets and OnSIDES. 
##############################################################################

lapply(c('Onsides','Opentargets'), function(dataset){
  
  datafile<-list.files(path='Data/', pattern='withgenetics', full.names=T)
  datafile=datafile[grepl(paste0(dataset), datafile)]
  withdrawn<-list.files(path='Data/', pattern='drugwarnings', full.names=T)
  withdrawn=withdrawn[grepl(paste0(dataset), withdrawn)]
  # Load SE-GPS constructed above
  Genescores_beta_genetic_CV=fread(paste0(scorefolder, 'All_genescoresum_across_all_predictors_',dataset,'_outcome_',outcome,'_sideeffectproject.txt.gz'), data.table=F)
  
  ADR_severity=fread(paste0('Data/Adr_severity_scores_phecodeX.txt'), data.table=F)
  dataset1<-fread(datafile, data.table=F)
  dataset1$ID=paste0(dataset1$drugname,'-',dataset1$gene,'-',dataset1$parentterm)
  
  #Format each drug dataset and severe side effects in different format. Opentargets
  #mapped to toxicity classes whereas OnSIDES at parentterm level
  
  if (dataset == 'Opentargets'){
    #the Open target warning dataset is mapped to toxicity classes - map these toxicity classes to phecode categories
    dataset1_withdrawn<-fread(withdrawn, data.table=F)
    colnames(dataset1_withdrawn)=gsub('warningType', 'warning_type',   colnames(dataset1_withdrawn))
    colnames(dataset1_withdrawn)=gsub('toxicityClass', 'warning_class',   colnames(dataset1_withdrawn))
    toxicity=c("gastrointestinal toxicity","cardiotoxicity","psychiatric toxicity","respiratory toxicity","hepatotoxicity","hematological toxicity",
               "metabolic toxicity","neurotoxicity","immune system toxicity","infectious disease","dermatological toxicity","musculoskeletal toxicity","nephrotoxicity","vascular toxicity","carcinogenicity")
    phecode_cat_tox= data.frame(warning_class = toxicity, category_warningclass = c("Gastrointestinal", "Cardiovascular", "Mental","Respiratory","Gastrointestinal","Blood/Immune","Endocrine/Metab", "Neurological","Blood/Immune","Infections","Dermatological","Muscloskeletal","Genitourinary","Cardiovascular","Neoplasms"))
    data_phecode_tox<-inner_join(phecode_cat_tox, dataset1_withdrawn,relationship = "many-to-many")
    dataset1_tox<-inner_join(dataset1, data_phecode_tox,relationship = "many-to-many")
    category=Genescores_beta_genetic_CV %>% distinct(category,category_full)
    dataset1_tox<-inner_join(dataset1_tox, category)
    dataset1_removese=dataset1 %>% select(-se,-senomi_phase4,-mi,-mi_phecode_phase4) %>% distinct()
  } else{
    dataset1_withdrawn<-fread(withdrawn, data.table=F)
    category=Genescores_beta_genetic_CV %>% distinct(category,category_full)
    # withdrawn drug-side effect dataset mapped to phecode parentterms
    dataset1_withdrawn1<-inner_join(dataset1_withdrawn,category, relationship = "many-to-many")
    dataset1_withdrawn2=dataset1_withdrawn1 %>% select(-pt_meddra_term) %>% mutate(se=1)  %>% mutate(category_warningclass=category_full)
    # combine withdrawn drug-side effect dataset with genetic predictors
    dataset1_tox<-inner_join(dataset1, dataset1_withdrawn2,relationship = "many-to-many")  
    dataset1_removese=dataset1 %>% select(-se,-senomi_phase4,-mi) %>% distinct()
    
  }
  # restrict to se toxicity phecodes present 
  lapply(c('senomi_phase4'), function(outcome){
    dataset1_tox_se=dataset1_tox %>% distinct(drugname, se,mi, senomi_phase4, parentterm, category_warningclass,category_full) %>% filter(.data[[outcome]]==1) %>% filter(category_warningclass==category_full)
    print(length(unique(dataset1_tox_se$parentterm)))
    print(length(unique(dataset1_tox$drugname)))
    # Format - Create a data frame with  all withdrawn drugs and side effect combinations
    dataset1_tox_se_dataformat <- expand.grid(drugname =unique(dataset1_tox_se$drugname), parentterm = unique(dataset1_tox_se$parentterm))
    cat_warninggroups<-dataset1_tox_se %>% distinct(parentterm, category_warningclass)
    # Combine with severe drug side effect data and side effect categories
    dataset1_tox_se_dataformat1<-left_join(dataset1_tox_se_dataformat, dataset1_tox_se[c('drugname','mi','se','senomi_phase4','parentterm')])
    dataset1_tox_se_dataformat2<-left_join(dataset1_tox_se_dataformat1, cat_warninggroups)
    dataset1_tox_se_dataformat2[is.na(dataset1_tox_se_dataformat2)]<-0
    colnames(dataset1)=gsub('category_full','category_warningclass',colnames(dataset1))
    #keep only se that are a drug warning
    dataset1_tox_se_dataformat2_gen<-inner_join(dataset1_removese, dataset1_tox_se_dataformat2,relationship = "many-to-many")
    print(nrow(dataset1_tox_se_dataformat2_gen %>% distinct(gene,drugname)) *length(unique(dataset1_tox_se_dataformat2_gen$parentterm))==dim(dataset1_tox_se_dataformat2_gen)[1])
    # Extract side SE-GPS and combine with severe drug dataset
      Genescores_beta_genetic_CV1= Genescores_beta_genetic_CV %>% distinct(drugname, gene, parentterm, genescoresum)  %>%arrange(genescoresum) %>% mutate(order=c(seq(1:length(genescoresum)))) %>% mutate(percent=order/length(genescoresum) *100) 
    Dataset_genescores1=inner_join(Genescores_beta_genetic_CV1, dataset1_tox_se_dataformat2_gen)
    write.table(Dataset_genescores1, paste0('Data/Withdrawndrugs_dataset_',dataset,'_outcome_',outcome,'.txt'), sep='\t', quote=F, row.names=F)
    
  })})



end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

