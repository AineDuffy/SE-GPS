# Side effect genetic-priority-score (SE-GPS)

<h2>Summary</h2>
<p>We created an in-silico side effect genetic priority score (SE-GPS) that can inform the likelihood of an adverse side effect across 15,139 genes and 499 phecode side effects. This score is constructed as a weighted sum using a multivariable mixed-effect regression model of four constructed genetic features with drug side effects to obtain the effect sizes from the association of each feature. These features include: 1) clinical variant evidence from ClinVar, HGMD and OMIM, consolidated into a single feature quantified as the number of overlapping entries; 2) single coding variants encompassing pLOF and missense single variants curated from Genebass and RAVAR; 3) Gene burden tests from Open Targets and RAVAR and 4) genome-wide association (GWA) loci, represented by two separate features: Locus2Gene and eQTL phenotype. This was applied to 19,422 protein-coding genes and 502 phecodes. We further incorporated the direction of genetic effect in a complementary score, the SE-GPS-DOE, to mimic the mechanism of the drug using predictions of loss-of-function (LOF) and gain-of-function (GOF) from LoGoFunc for clinical variants and QTL estimates for GWAS phenotypes. Positive directional scores reflect the likelihood of an adverse side effect from target inhibition, and a negative SE-GPS-DOE reflects target activation.</p>

<h2>Study</h2>
<p>For further details about the methods and analysis behind this study, please see our paper: Duffy, A et al. Development of a Genetic Priority Score to Predict Drug Side Effects Using Human Genetic Evidence. Submitted.</p>
These scores can be explored further at: https://rstudio-connect.hpc.mssm.edu/sideeffect-geneticpriorityscore/

<h2>Instructions</h2>
This repository contains all the necessary scripts to create the SE-GPS and the SE-GPS-DOE from the final Open Target and OnSIDES genetic datasets, as well as to reproduce all figures from the manuscript. All scripts should take less than 20mins. The data can be accessed at https://zenodo.org/uploads/15334136. We provide three scripts, which should be run in the following order:

<b></b>
1. Create_scores.R  - Creates the training and test Open Target datasets, runs the mixed model, and generates the SE-GPS across Open Targets, OnSIDES, and all genes.
2. Create_scores_DOE.R - Creates the training and test Open Target datasets incorporating direction of effect, runs the mixed model, and generates the SE-GPS-DOE across Open Targets, OnSIDES, and all genes.
3. Analysis_Figures.R  -  Code to run all analyses and figures  

This project uses renv to manage the R environment and ensure that all dependencies are properly handled. After cloning the repository, run the following command in your R session to restore the project’s environment from the renv.lock file: renv::restore().
