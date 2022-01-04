# Current array
# array.id = as.integer(commandArgs()[length(commandArgs())])

#' This is revised version of the simulation code
#' This code will be used for simulating the type I error
#' 
# install.packages("MASS")
# install.packages("gdata")
# install.packages("CompQuadForm")
# install.packages("dglm")
# install.packages("lubridate")

library(gridExtra)
library(reshape2)
library(lubridate)
library(tidyverse)

setwd("/n/holystore01/LABS/xlin/Lab/xihaoli/DOHG/")

###### Predefined global parameter
#' Number of simulations, number of subjects, Number of phenotypes and MAF
#' (independent of covariates)

Nsub_array   <- c(5000, 10000)
Npheno_array <- c(4, 6, 8, 10, 12)
maf_array    <- c(0.25)
gamma_array  <- c(NA, 0.15, 0.55, 0.95)
job_slice_array <- seq(1, 100)
Nsim         <- 10000

param_comb <- expand.grid(Nsub   = Nsub_array,
                          Npheno = Npheno_array,
                          maf    = maf_array,
                          gamma  = gamma_array,
                          JobSlice = job_slice_array)

results <- matrix(NA, 40, 9)
colnames(results) <- c("Nsub", "Npheno", "maf", "gamma", "JobSlice", "SMAT", "SMAT_X2", "CCT", "CCT_X2")
for (array.id in 1:40) {
  Nsub   <- param_comb[array.id, "Nsub"]
  Npheno <- param_comb[array.id, "Npheno"]
  maf    <- param_comb[array.id, "maf"]
  gamma  <- param_comb[array.id, "gamma"]
  JobSlice <- param_comb[array.id, "JobSlice"]
  
  results[array.id, 1] <- Nsub
  results[array.id, 2] <- Npheno
  results[array.id, 3] <- maf
  results[array.id, 4] <- gamma
  results[array.id, 5] <- JobSlice
  
  file_SMAT_pvalue               <- paste("pvalue", "SMAT", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_includeXsquar_SMAT_pvalue <- paste("pvalue", "includeXsquare_SMAT", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_CCT_pvalue                <- paste("pvalue", "CCT", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_includeXsquare_CCT_pvalue <- paste("pvalue", "includeXsquare_CCT", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  
  file_SMAT_time                <- paste("time", "SMAT_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_includeXsquare_SMAT_time <- paste("time", "includeXsquare_SMAT_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_CCT_time                 <- paste("time", "CCT_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_includeXsquare_CCT_time  <- paste("time", "includeXsquare_CCT_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_simulation_time          <- paste("time", "simulation_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  
  if (file.exists(paste0("Result/New_4_6_Simulation_ExcludeGAMUT/", file_SMAT_pvalue, ".RData")) & 
      file.exists(paste0("Result/New_4_6_Simulation_ExcludeGAMUT/", file_includeXsquar_SMAT_pvalue, ".RData")) & 
      file.exists(paste0("Result/New_4_6_Simulation_ExcludeGAMUT/", file_CCT_pvalue, ".RData")) & 
      file.exists(paste0("Result/New_4_6_Simulation_ExcludeGAMUT/", file_includeXsquare_CCT_pvalue, ".RData"))) {
    
    load(paste0("Result/New_4_6_Simulation_ExcludeGAMUT/", file_SMAT_pvalue, ".RData"))
    load(paste0("Result/New_4_6_Simulation_ExcludeGAMUT/", file_includeXsquar_SMAT_pvalue, ".RData"))
    load(paste0("Result/New_4_6_Simulation_ExcludeGAMUT/", file_CCT_pvalue, ".RData"))
    load(paste0("Result/New_4_6_Simulation_ExcludeGAMUT/", file_includeXsquare_CCT_pvalue, ".RData"))
    
    results_temp <- data.frame(pval_SMAT_adj_dglm, pval_SMAT_includeXsquare_adj_dglm, pval_CCT_adj_dglm, pval_includeXsquare_CCT_adj_dglm)
    results[array.id, c(6:9)] <- apply(results_temp, 2, function(x) {sum(x < 0.01) / dim(results_temp)[1]})
  }
}

results <- data.frame(results)
