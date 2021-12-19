library(tidyverse)
install.packages("qqman")
library(qqman)

setwd("~/Dropbox/Emory Courses/DOHG")
source("Analysis/result_analysis_source_code.R")

load("Result/New_4_6_Simulation/pvalue_CCT_adj_dglm_10000_subj_10000_pheno_4_MAF_0.15.RData")
pval_CCT <- pval_CCT_adj_dglm

load("Result/New_4_6_Simulation/pvalue_includeXsquare_CCT_adj_dglm_10000_subj_10000_pheno_4_MAF_0.15.RData")
pval_includeXsquare_CCT <- pval_includeXsquare_CCT_adj_dglm

load("Result/New_4_6_Simulation/pvalue_SMAT_adj_dglm_10000_subj_10000_pheno_4_MAF_0.15.RData")
pval_SMAT <- pval_SMAT_adj_dglm

###
hist(pval_CCT,
     main = "Type 1 Error of Cauchy Combination Test \n4 Phenotypes on 10,000 subjects \n10,000 Simulations")
qqunif.plot(as.numeric(pval_CCT), 
            main = "Type 1 Error of Cauchy Combination Test \n4 Phenotypes on 10,000 subjects \n10,000 Simulations With 95% Confidence Interval",
            col = "#012169")

###
hist(pval_includeXsquare_CCT,
     main = "Type 1 Error of Cauchy Combination Test with X^2 \n4 Phenotypes on 10,000 subjects \n10,000 Simulations")
qqunif.plot(as.numeric(pval_includeXsquare_CCT), 
            main = "Type 1 Error of Cauchy Combination Test with X^2 \n4 Phenotypes on 10,000 subjects \n10,000 Simulations With 95% Confidence Interval",
            col = "#012169")

###
hist(pval_SMAT,
     main = "Type 1 Error of SMAT \n4 Phenotypes on 10,000 subjects \n10,000 Simulations")
qqunif.plot(as.numeric(pval_SMAT), 
            main = "Type 1 Error of SMAT \n4 Phenotypes on 10,000 subjects \n10,000 Simulations With 95% Confidence Interval",
            col = "#012169")
