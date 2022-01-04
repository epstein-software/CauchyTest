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

Nsub_array   <- c(5000)
# Nsub_array   <- c(5000, 10000)
# Npheno_array <- c(4, 6, 8, 10, 12)
Npheno_array <- c(4, 6, 8)
maf_array    <- c(0.25)
gamma_array  <- c(NA, 0.15, 0.55, 0.95)
job_slice_array <- seq(1, 100)
Nsim         <- 10000

param_comb <- expand.grid(Nsub   = Nsub_array,
                          Npheno = Npheno_array,
                          maf    = maf_array,
                          gamma  = gamma_array,
                          JobSlice = job_slice_array)

############ ---------- Construct the Overall Dataframe ----------- ############

full_df <- data.frame()
for (array.id in 1:nrow(param_comb)) {
  
  Nsub   <- param_comb[array.id, "Nsub"]
  Npheno <- param_comb[array.id, "Npheno"]
  maf    <- param_comb[array.id, "maf"]
  gamma  <- param_comb[array.id, "gamma"]
  JobSlice <- param_comb[array.id, "JobSlice"]
  
  file_SMAT_pvalue               <- paste("pvalue", "SMAT", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_includeXsquar_SMAT_pvalue <- paste("pvalue", "includeXsquare_SMAT", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_CCT_pvalue                <- paste("pvalue", "CCT", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_includeXsquare_CCT_pvalue <- paste("pvalue", "includeXsquare_CCT", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  
  file_SMAT_time                <- paste("time", "SMAT_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_includeXsquare_SMAT_time <- paste("time", "includeXsquare_SMAT_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_CCT_time                 <- paste("time", "CCT_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_includeXsquare_CCT_time  <- paste("time", "includeXsquare_CCT_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  file_simulation_time          <- paste("time", "simulation_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
  
  # assign(file_SMAT_time,
  #        get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_SMAT_time, ".RData", sep = ""))))
  # assign(file_includeXsquare_SMAT_time,
  #        get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_includeXsquare_SMAT_time, ".RData", sep = ""))))
  # assign(file_CCT_time,
  #        get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_CCT_time, ".RData", sep = ""))))
  # assign(file_includeXsquare_CCT_time,
  #        get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_includeXsquare_CCT_time, ".RData", sep = ""))))
  # assign(file_simulation_time,
  #        get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_simulation_time, ".RData", sep = ""))))
  # 
  # assign(file_SMAT_pvalue,
  #        get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_SMAT_pvalue, ".RData", sep = ""))))
  # assign(file_includeXsquar_SMAT_pvalue,
  #        get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_includeXsquar_SMAT_pvalue, ".RData", sep = ""))))
  # assign(file_CCT_pvalue,
  #        get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_CCT_pvalue, ".RData", sep = ""))))
  # assign(file_includeXsquare_CCT_pvalue,
  #        get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_includeXsquare_CCT_pvalue, ".RData", sep = ""))))
  
  temp_DF <- data.frame(
    "Nsub"          = Nsub,
    "Npheno"        = Npheno,
    "maf"           = maf,
    "gamma"         = gamma,
    "JobSlice"            = JobSlice,
    "SMAT"                = get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_SMAT_pvalue, ".RData", sep = ""))),
    "SMAT_includeXsquare" = get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_includeXsquar_SMAT_pvalue, ".RData", sep = ""))),
    "CCT"                 = get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_CCT_pvalue, ".RData", sep = ""))),
    "CCT_includeXsquare"  = get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_includeXsquare_CCT_pvalue, ".RData", sep = "")))
    #    "SMAT_Time"             = get(load(paste("Result/New_4_6_Power_Simulation/", file_SMAT_time, ".RData", sep = ""))),
    #    "CCT_Time"              = get(load(paste("Result/New_4_6_Power_Simulation/", file_CCT_time, ".RData", sep = ""))),
    #    "CCT_includeXsquare_Time" = get(load(paste("Result/New_4_6_Power_Simulation/", file_includeXsquare_CCT_time, ".RData", sep = ""))),
    #    "Overall_Time" = get(load(paste("Result/New_4_6_Power_Simulation/", file_simulation_time, ".RData", sep = "")))
  )
  
  full_df <- rbind(full_df, temp_DF)
}