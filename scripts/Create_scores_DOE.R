# Load required libraries
library(data.table)
library(dplyr)

set.seed(81)

start.time <- Sys.time()
print(start.time)

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
#######################################
# 1. Create 5 fold CV training sets 
#######################################

#load the Open target drug-genetic dataset with directional predictors. 
Opentarget_doe=fread(paste0('Data/','/Opentargetsgeneticdataset_filtered_5_All_DOE_matchmechanism.txt.gz'), data.table=F)
Opentarget_doe$gene_phenotype=paste0(Opentarget_doe$gene,'-',Opentarget_doe$parentterm)

# Create five equally sized groups of randomly assigned unique gene-phenotype pairs -  5 cross validation 20% sets
# Each group is test set and training set will be remaining 80%
uniquegene_phenotype=Opentarget_doe%>% distinct(gene_phenotype)
uniquegene_phenotype_ss <- sample(rep(1:5, diff(floor(nrow(uniquegene_phenotype) * c(0, 0.2, 0.4, 0.6, 0.8, 1)))))
splitdata_gp <- split(uniquegene_phenotype,uniquegene_phenotype_ss)

for(i in 1:length(splitdata_gp)){
  data_split=splitdata_gp[[i]]
  test_set <- Opentarget_doe[Opentarget_doe$gene_phenotype%in% data_split$gene_phenotype, ]
  filename=paste0('results/construct_scores/',dataset,'_geneticdataset_CVsample',i,'_DOE_matchmechanism.txt.gz')
  print(filename)
  write.table(test_set, gzfile(filename), sep='\t', quote=F, row.names=F)
}

# unload and reload required libraries
unloadNamespace("data.table")
unloadNamespace("dplyr")

# Load required libraries
library("Matrix")
library("Rcpp")
library("StanHeaders")
library("rstan")
options(mc.cores = parallel::detectCores())
library("lme4")
library("brms")
library(data.table)
library(parallel)

#######################################
# 2. Run mixed model in 80% training set 
#######################################

dataset='Opentargets';independent='onside';outcome='senomi_phase4'
geneticpredictors= c('clinicalvariant', 'gwastrait','geneburden','singlevar')
Predictor_names= data.frame(from = geneticpredictors, to = c("clinicalvariant_DOE","gwastrait_doe","geneburden_doe","singlevar_doe"))

#load the Open target drug-genetic dataset with directional predictors
#Directional dataset restricted to drugs with either an activator or inhibitor mechanism 
doedatafile<-list.files(path='Data/', pattern='DOE', full.names=T)
doedatafile=doedatafile[grepl(paste0('Opentargetsgeneticdataset_filtered'), doedatafile)]

doe_dataset<-fread(doedatafile, data.table=F)
doe_dataset1<-doe_dataset[c('drugname','gene','parentterm','se','senomi_phase4','moa')]
doe_dataset1_full<-merge(doe_dataset1, opentarget_fulldataset)
print(head(doe_dataset1_full))

# Combine the mixed effects results model across 5 cross-validation (CV) samples
Mixedmodel_all<-do.call(rbind, lapply(c(paste0('CVsample',rep(1:5))), function(CVsample){
  CVdataset=fread(paste0(scorefolder,dataset,'_geneticdataset_',CVsample,'_DOE_matchmechanism.txt.gz'),data.table=F) #20% training Open target dataset
  CVdataset$ID=paste0(CVdataset$drugname,'-',CVdataset$gene,'-',CVdataset$parentterm)
  dataset_sampled80=subset(doe_dataset1_full, !(ID %in% CVdataset$ID))#80% training drug dataset
  geneticpredictors= c('clinicalvariant', 'gwastrait','geneburden','singlevar')
  ADR_severity1=fread(paste0('Data/Adr_severity_scores_phecodeX.txt'), data.table=F)
  dataset_sampled80_ADR=merge(dataset_sampled80, ADR_severity1, all.x=T)
  dataset_sampled80_ADR$genepheno=paste0(dataset_sampled80_ADR$gene,'-', dataset_sampled80_ADR$parentterm)
  dataset_sampled80_ADR$severity_score_round=round(dataset_sampled80_ADR$severity_score*100)
  Predictors=paste0(paste(geneticpredictors,collapse='+'), '+category + moa ') # include mechanism of action (moa) as covariate as well
  ##Run mixedmodel regression across each CV training dataset to get weights
  mod_se_mixed <- glmer(as.formula(paste0(outcome, " ~ ", paste(Predictors), " +(1 | drugname)")),
                        family = "binomial",
                        data = dataset_sampled80_ADR,
                        control = glmerControl(optimizer = "bobyqa"),nAGQ=0L,
                        weights=severity_score_round)
  
  mod_se_mixed_results1 <- cbind.data.frame(beta = summary(mod_se_mixed)$coefficient[, c(1)],
                                            lowerCI = summary(mod_se_mixed)$coefficient[, c(1)] - 1.96 * summary(mod_se_mixed)$coefficient[, c(2)],
                                            upperCI = summary(mod_se_mixed)$coefficient[, c(1)] + 1.96 * summary(mod_se_mixed)$coefficient[, c(2)],
                                            P.val = summary(mod_se_mixed)$coefficient[, c(4)])
  
  mod_se_mixed_results1$Predictor=rownames(mod_se_mixed_results1)
  mod_se_mixed_results1$randomeffect='drug'
  mod_se_mixed_results1$CV=CVsample
  
  write.table(mod_se_mixed_results1, paste0(scorefolder, '/Mixedeffect_weighted_model_ADR_severity_',dataset,'_outcome_',outcome,'_',CVsample,'_DOE_matchmechanism.txt'), sep='\t',quote=F)
  
  return(mod_se_mixed_results1)
  
}))

print(paste('load all packges...'))
print(paste(head(Mixedmodel_all)))
library(plyr)
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(DataCombine)
library(logistf)
library(parallel)


Sup_table_weights<-Mixedmodel_all %>% dplyr::mutate(CI=paste0(round(lowerCI,2), ' - ', round(as.numeric(upperCI),2)), beta=round(beta,2)) %>%
  select(CV, Predictor, beta, CI, P.val)  %>% dplyr::rename(`P-value`=`P.val`) %>% arrange(CV)

write.table(Sup_table_weights, paste0(tables_dir,'/Table_S3_Mixedmodel_',dataset,'_multivariatemodel_results.txt') , sep='\t', quote=F, row.names=F)



##############################################################################
# 3. Create DOE scores - apply betas from mixed model as weights to CV test set
##############################################################################

Genescore_sum<-lapply(c(paste0('CVsample',rep(1:5))), function(CVsample){
  
  mixedmodelresults=fread(paste0(scorefolder, '/Mixedeffect_weighted_model_ADR_severity_',dataset,'_outcome_',outcome,'_',CVsample,'_DOE_matchmechanism.txt'), data.table=F)
  mixedmodelresults=mixedmodelresults[!grepl('category|Intercept|moainhibitor', mixedmodelresults$Predictor),]
  mixedmodelresults <- FindReplace(data = mixedmodelresults, Var = "Predictor", replaceData = Predictor_names,
                                   from = "from", to = "to", exact = TRUE)
  
  # Use weights from Firth and create score using remaining 20% of data for each CV
  Firth_weights=as.data.frame(t(as.numeric(mixedmodelresults$beta)))
  names(Firth_weights)=mixedmodelresults$Predictor
  #extract gene, parentterm and genetic predictor columns
  CVdataset=fread(paste0(scorefolder,dataset,'_geneticdataset_',CVsample,'_DOE_matchmechanism.txt.gz'),data.table=F) #20% training Open target dataset
  CVdataset_gene_phenotype=CVdataset[, grepl(paste0('\\bgene\\b|parentterm|moa|',paste0(Predictor_names$to, collapse='\\b|')), colnames(CVdataset))] %>% distinct() 
  
  # Multiply each predictor value in the test set by its corresponding beta coefficient weight
  Genescores_beta_genetic_values=mapply("*", CVdataset_gene_phenotype[intersect(names(CVdataset_gene_phenotype), names(Firth_weights))],
                                        Firth_weights[intersect(names(CVdataset_gene_phenotype), names(Firth_weights))])
  gene_parentterm=CVdataset_gene_phenotype %>% select(gene, parentterm, moa)
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #Construct SE-GPS-DOE as the sum across all predictor values for each gene-phenotype pair.
  # Negative scores from GOF evidence reflect target activation, positive scores from LOF evidence reflect target inhibition
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm","moa")])
  write.table(Genescores_beta_genetic, paste0(scorefolder, 'Genescore_sum_across_predictor_',dataset,'_outcome_',outcome,'_',CVsample,'_DOE.txt'), sep='\t', row.names=F, quote=F )
  
})

#Choose the CV sample with the max OR for validation set

max_filetype<-do.call(rbind,lapply(c(paste0('CVsample',seq(1:5))), function(CVsample){   
  CVdataset=fread(paste0(scorefolder,dataset,'_geneticdataset_',CVsample,'_DOE_matchmechanism.txt.gz'),data.table=F) #20% training Open target dataset
  genescorefile=fread(paste0(scorefolder, 'Genescore_sum_across_predictor_',dataset,'_outcome_',outcome,'_',CVsample,'_DOE.txt'), data.table=F)
  # combine genescore sum file with drug data for each test set 
  Dataset_genescores=inner_join(CVdataset[c('drugname','gene','parentterm','moa','category','category_full','mi','se','senomi_phase4')] ,genescorefile, by=c('gene','parentterm') ,relationship = "many-to-many")
  #run logistic model 
  model1 =glm(as.formula(paste0(outcome, ' ~ genescoresum + category')), data=Dataset_genescores,family = 'binomial')
  mod_output<-rbind(cbind.data.frame(CV=CVsample,outcome=outcome,OR=exp(summary(model1)$coefficient[2,1]),lowerCI=exp(summary(model1)$coefficient[2,1]-(1.96* summary(model1)$coefficient[2,2])),upperCI=exp(summary(model1)$coefficient[2,1]+(1.96* summary(model1)$coefficient[2,2])),P.val=summary(model1)$coefficient[2,4]))
  return(mod_output)
}))

## take CV with the max OR

Sup_table<-max_filetype %>% mutate(CI=paste0(round(lowerCI,2), ' - ', round(as.numeric(upperCI),2)), OR=round(OR,2)) %>%
  select(CV, OR, CI, P.val)  %>% dplyr::rename(`P-value`=`P.val`) %>% arrange(desc(OR))

write.table(Sup_table, paste0(tables_dir, 'Table_S4_associations_across_CVtest_sets_DOE.txt'), sep='\t', quote=F, row.names=F)

## For each CV-fold, merge the SE-GPS-DOE with the corresponding CV 20% drug dataset at the gene-phenotype level to combine the SE-GPS-DOE with drug side effect information
## combine into a single data frame (100% of data)

Combine_genescores<-lapply(paste0('CVsample', seq(1:5)), function(CVsample){
  CVdataset=fread(paste0(scorefolder,dataset,'_geneticdataset_',CVsample,'_DOE_matchmechanism.txt.gz'),data.table=F) #20% training Open target dataset
  genescorefile=fread(paste0(scorefolder, 'Genescore_sum_across_predictor_',dataset,'_outcome_',outcome,'_',CVsample,'_DOE.txt'), data.table=F)
  Genescore_sumfile_20_drugs<-inner_join(CVdataset[c('drugname','gene','parentterm','moa','category','category_full','mi','se','senomi_phase4')], genescorefile ) %>% distinct()
  Genescore_sumfile_20_drugs$CV=CVsample
  return(Genescore_sumfile_20_drugs)
})
Combine_genescores1<-do.call(rbind.fill,Combine_genescores)
write.table(Combine_genescores1, gzfile(paste0(scorefolder, 'All_genescoresum_across_all_predictors_',dataset,'_outcome_',outcome,'_sideeffectproject_DOE.txt.gz')), sep='\t',quote=F,row.names=F)


##############################################################################
# 4. Apply method to OnSIDES dataset 
##############################################################################

lapply(c('Onsides'), function(valdataset){
  #Calculate SE-GPS-DOE in OnSIDES dataset using the OT model weights that produced the highest OR 
  OT_CV=fread(paste0(tables_dir, 'Table_S4_associations_across_CVtest_sets_DOE.txt'), data.table=F )
  samplecv=head(OT_CV[1],1)[[1]]
  # Load mixed DOE model coefficients (predictor weights) from Open Targets for the selected CV fold
  mixedmodelresults=fread(paste0(scorefolder, '/Mixedeffect_weighted_model_ADR_severity_',dataset,'_outcome_',outcome,'_',samplecv,'_DOE_matchmechanism.txt'), data.table=F)
  mixedmodelresults=mixedmodelresults[!grepl('category|Intercept|moainhibitor', mixedmodelresults$Predictor),]
  mixedmodelresults <- FindReplace(data = mixedmodelresults, Var = "Predictor", replaceData = Predictor_names,
                                   from = "from", to = "to", exact = TRUE)
  Firth=mixedmodelresults %>% select(Predictor, beta) 
  valdataset_dir<-list.files(path='Data', pattern='Onsidesgeneticdataset', full.names=T)
  # Load the OnSIDES directional validation dataset
  valdataset_dir=valdataset_dir[grepl(paste0('_DOE'), valdataset_dir)]
  Validation_dataset=fread(paste0(valdataset_dir),data.table=F) 
  # Extract beta coefficients 
  Firth_weights=as.data.frame(t(Firth$beta))
  names(Firth_weights)=Firth$Predictor
  # Select unique gene-phenotype predictor combinations from the validation dataset
  Validation_dataset_gene_phenotype=Validation_dataset[, grepl(paste0('\\bgene\\b|parentterm|moa|',paste0(Predictor_names$to, collapse='\\b|')), colnames(Validation_dataset))] %>% distinct() 
  # Multiply each predictor value in the validation set by its corresponding OT beta coefficient weight
  gene_parentterm=Validation_dataset %>% distinct(gene, parentterm, moa)
  Genescores_beta_genetic_values=mapply("*", Validation_dataset_gene_phenotype[intersect(names(Validation_dataset_gene_phenotype), names(Firth_weights))],
                                        Firth_weights[intersect(names(Validation_dataset_gene_phenotype), names(Firth_weights))])
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #Construct SE-GPS-DOE as the sum across all predictor values for each gene-phenotype pair
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm","moa")])
  Genescores_beta_genetic_drugs=inner_join(Validation_dataset[c('drugname','gene','parentterm','moa','category','category_full','mi','se','senomi_phase4')],Genescores_beta_genetic, by=c('gene', 'parentterm','moa'),relationship = "many-to-many")
  write.table(Genescores_beta_genetic_drugs, paste0(scorefolder, 'All_genescoresum_across_all_predictors_',valdataset,'_outcome_',outcome,'_sideeffectproject_DOE.txt.gz'), sep='\t', row.names=F, quote=F )
  Genescores_beta_genetic_genept=Genescores_beta_genetic %>% distinct(gene, parentterm, genescoresum, moa) 
  genescorefile_drugs=inner_join(Validation_dataset[c('drugname','gene','parentterm','moa','category','mi','se','senomi_phase4')],Genescores_beta_genetic_genept, by=c('gene', 'parentterm'),relationship = "many-to-many")
  write.table(genescorefile_drugs, paste0(scorefolder, 'All_genescoresum_',valdataset,'_outcome_',outcome,'_sideeffectproject_DOE.txt.gz'), sep='\t', row.names=F, quote=F )
  
})

##############################################################################
# 5. Apply method to All genes  
##############################################################################

lapply(c('Allgenes'), function(valdataset){
  
  # Calculate SE-GPS-DOE in all genes dataset (19000 +genes) using the OT mixedmodel weights which gave the max OR 
  OT_CV=fread(paste0(tables_dir, 'Table_S4_associations_across_CVtest_sets_DOE.txt'), data.table=F )
  samplecv=head(OT_CV[1],1)[[1]]
  # Load mixed model coefficients (predictor weights) from Open Targets for the selected CV fold
  mixedmodelresults=fread(paste0(scorefolder, '/Mixedeffect_weighted_model_ADR_severity_',dataset,'_outcome_',outcome,'_',samplecv,'_DOE_matchmechanism.txt'), data.table=F)
  # Remove intercept and covariate coefficients
  mixedmodelresults=mixedmodelresults[!grepl('category|Intercept|moainhibitor', mixedmodelresults$Predictor),]
  mixedmodelresults <- FindReplace(data = mixedmodelresults, Var = "Predictor", replaceData = Predictor_names,
                                   from = "from", to = "to", exact = TRUE)
  # Extract beta coefficients 
  Firth=mixedmodelresults %>% select(Predictor, beta) 
  Firth_weights=as.data.frame(t(Firth$beta))
  names(Firth_weights)=Firth$Predictor
  # Load dataset with all gene-phenotype pairs and predictor values
  phenotype_allgene<-fread(paste0('Data/Allgenes_dataset_phenotype_predictors_stringentfilters_sideeffect_project.txt.gz'), data.table=F)
  
  #Gain of function predictions
  DOEGOF=fread(paste0('Data/Predictors_with_direction_effect_sideeffect_project_gof_predvalue.txt'),data.table=F)
  #Collapse from phecode to parentterm level
  DOEGOF1=as.data.frame(DOEGOF %>% mutate(parentterm=as.character(trunc(as.numeric(Phecode)))) %>%
                          group_by( gene, parentterm,category_full) %>% 
                          dplyr::summarise(omim_predicted:=max(omim_predicted),eqtl_Pheno_DOE:=max(eqtl_Pheno_DOE) ,hgmd_predicted:=max(hgmd_predicted),
                                           genebass_predicted:=max(genebass_predicted), OTgeneburden_doe:=max(OTgeneburden_doe), OTG_DOE:=max(OTG_DOE),
                                           clinvar_predicted:=max(clinvar_predicted), ravar_genburden_DOE:=max(ravar_genburden_DOE), ravar_snp_predicted:=max(ravar_snp_predicted) ) %>% distinct())
  DOEGOF1$moa='activator'
  DOEGOF1[, 4:12] <- lapply(DOEGOF1[, 4:12], function(x) ifelse(x !=0, 'GOF', x))
  #Loss of function predictons
  DOELOF=fread(paste0('Data/Predictors_with_direction_effect_sideeffect_project_lof_predvalue.txt'),data.table=F)
  #Collapse from phecode to parentterm level
  DOELOF1=as.data.frame(DOELOF %>% mutate(parentterm=as.character(trunc(as.numeric(Phecode)))) %>% group_by( gene, parentterm,category_full) 
                        %>% dplyr::summarise(omim_predicted:=max(omim_predicted),eqtl_Pheno_DOE:=max(eqtl_Pheno_DOE),hgmd_predicted:=max(hgmd_predicted),genebass_predicted:=max(genebass_predicted), OTgeneburden_doe:=max(OTgeneburden_doe), 
                                             OTG_DOE:=max(OTG_DOE),clinvar_predicted:=max(clinvar_predicted), ravar_genburden_DOE:=max(ravar_genburden_DOE), ravar_snp_predicted:=max(ravar_snp_predicted) ) %>% distinct())
  
  DOELOF1$moa='inhibitor'
  DOELOF1[, 4:12] <- lapply(DOELOF1[, 4:12], function(x) ifelse(x !=0, 'LOF', x))
  # Combine GOF and LOF predictions, convert GOF = -1, LOF = 1, NA = 0
  doe_both<-rbind(DOELOF1,DOEGOF1)
  doe_both$parentterm=as.character(doe_both$parentterm)
  doe_both[doe_both == "LOF"] <- 1
  doe_both[doe_both == "GOF"] <- -1
  doe_both[is.na(doe_both)] <- 0
  
  doe_predictors=c('clinvar_predicted','omim_predicted','hgmd_predicted','OTgeneburden_doe','ravar_genburden_DOE','ravar_snp_predicted','genebass_predicted','eqtl_Pheno_DOE','OTG_DOE')
  doe_both[,doe_predictors]<-lapply(doe_both[,doe_predictors], as.numeric)
  clincicalpred_DOE=c('clinvar_predicted','omim_predicted','hgmd_predicted')
  ## Make composite clinical Variant predictor 
  doe_both$clinicalvariant_DOE = rowSums(doe_both[, clincicalpred_DOE])
  
  # Create binary predictors for geneburden, single variant and gwas trait
  geneburden_doe=c('OTgeneburden_doe','ravar_genburden_DOE')
  doe_both$geneburden_doe = ifelse(rowSums(doe_both[, geneburden_doe]) != 0, 1, 0)

  singlevar_doe=c('ravar_snp_predicted','genebass_predicted')
  doe_both$singlevar_doe = sign(rowSums(doe_both[, singlevar_doe]))

  gwastrait_doe=c('OTG_DOE','eqtl_Pheno_DOE')
  doe_both$gwastrait_doe = sign(rowSums(doe_both[, gwastrait_doe]))
  
  Validation_dataset_gene_phenotype=doe_both[, grepl(paste0('\\bgene\\b|parentterm|moa|',paste0(Predictor_names$to, collapse='\\b|')), colnames(doe_both))] %>% distinct() 
  
  # Multiply each predictor value in the All genes data set by its corresponding OT beta coefficient weight
  gene_parentterm=Validation_dataset_gene_phenotype %>% distinct(gene, parentterm, moa)
  Genescores_beta_genetic_values=mapply("*", Validation_dataset_gene_phenotype[intersect(names(Validation_dataset_gene_phenotype), names(Firth_weights))],
                                        Firth_weights[intersect(names(Validation_dataset_gene_phenotype), names(Firth_weights))])
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #Construct SE-GPS-DOE as the sum across all predictor values for each gene-phenotype pair
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm","moa")])
  # include phecode categories
  phecodef=read.csv('Data/phecodeX_info.csv')
  phecodef = phecodef %>% dplyr::rename(category_full=category) %>% mutate(category=gsub('_.*', '',phecode), phecode=gsub('.*_', '',phecode)) %>%
    mutate(parentterm=as.character(trunc(as.numeric(phecode)))) %>%  select(parentterm, category_full, category ) %>% distinct() 
  Genescores_beta_genetic_drugs=inner_join(phecodef,Genescores_beta_genetic, by=c( 'parentterm'),relationship = "many-to-many")
  write.table(Genescores_beta_genetic_drugs, paste0(scorefolder, 'All_genescoresum_across_all_predictors_',valdataset,'_outcome_',outcome,'_DOE.txt.gz'), sep='\t', row.names=F, quote=F )
})

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)


