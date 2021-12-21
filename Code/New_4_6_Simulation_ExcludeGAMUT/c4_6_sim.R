# Current array
array.id = as.integer(commandArgs()[length(commandArgs())])

#' This is revised version of the simulation code
#' This code will be used for simulating the type I error
#' 
# install.packages("MASS")
# install.packages("gdata")
# install.packages("CompQuadForm")
# install.packages("dglm")
# install.packages("lubridate")

library(MASS)
library(gdata)
library(CompQuadForm)
library(dglm)
library(lubridate)

setwd("/mnt/EpsteinFSS/data/sbian/CauchyTest/")
# setwd("~/Dropbox/Emory Courses/DOHG/")

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

paste("----------Start:", array.id, "out of", dim(param_comb)[1], "----------", sep = " ")

Nsub   <- param_comb[array.id, "Nsub"]
Npheno <- param_comb[array.id, "Npheno"]
maf    <- param_comb[array.id, "maf"]
gamma  <- param_comb[array.id, "gamma"]
JobSlice <- param_comb[array.id, "JobSlice"]

source('Code/New_4_6_Simulation_ExcludeGAMUT/c4_6_sim_source_code.R') 

#' Seed used for simulating the coefficients, the coefficients will be fixed in
#' the beginning. Instead of random simulating like this, we could 
int_array            <- array(NA,dim=c(Npheno))
main_beta_g_array    <- array(NA,dim=c(Npheno))
main_gamma_con_array <- array(NA,dim=c(Npheno))

set.seed( (array.id-1)%%40 + 1)

for (sim_Npheno in 1:Npheno) {
  
  # randomly-generated intercept for phenotypes
  int_array[sim_Npheno] <- rnorm(1, 0, 5) 
  
  # randomly-generated main effect of SNP on phenotypes
  main_beta_g_array[sim_Npheno] <- runif(1, 0, 0.2)
  
  # randomly-generate main effect of continuous covariate on genotype
  # Random simulate gamma is gamma == NA
  if (is.na(gamma)) {
    main_gamma_con_array[sim_Npheno] <- runif(1, 0, 0.2)
  } else {
    main_gamma_con_array[sim_Npheno] <- gamma
  }
}

#' Number of off-diagonal elements in covariance matrix of phenotype
#' (within subject)
Num_off_diag <- gamma(Npheno+1)/(gamma(3)*gamma(Npheno+1-2)) 

# Initialize arrays/vectors to contain p-values and time 
# pval_GAMuT_ProjMatrix_adj_dglm <- array(NA,dim=c(Nsim))
# pval_GAMuT_LineKernel_adj_dglm <- array(NA,dim=c(Nsim))
pval_SMAT_adj_dglm                <- array(NA,dim=c(Nsim))
pval_SMAT_includeXsquare_adj_dglm <- array(NA,dim=c(Nsim))
pval_CCT_adj_dglm                 <- array(NA,dim=c(Nsim))
pval_includeXsquare_CCT_adj_dglm  <- array(NA,dim=c(Nsim))

# GAMuT_ProjMatrix_time = c()
# GAMuT_LineKernel_time = c()
SMAT_time = c()
SMAT_includeXsquare_time = c()
CCT_time = c()
CCT_includeXsquare_time = c()
simulation_time = c()

for(isim in 1:Nsim) {
  
  set.seed(isim + 10000*(JobSlice - 1))
  
  start_time <- Sys.time()
  
  #' Each subject has their own X and W. Within each subject, X and W are fixed,
  #' but the coefficients vary again the phenotypes. The phenotypes vary by 
  #' randomly simulating the error term. This step results in (Nsub by Npheno)
  #'  simulated Y 
  X     <- rbinom(Nsub, 2, maf)
  W_con <- rnorm(Nsub, 0, 1)

  sim_Y     <- matrix(, nrow = Nsub, ncol = Npheno)
  sim_error <- matrix(, nrow = Nsub, ncol = Npheno)
  
  for (ipheno in 1:Npheno) {
    
    e_error <- rnorm(Nsub, 0, 1)
    Y <- int_array[ipheno] + main_beta_g_array[ipheno]*X + main_gamma_con_array[ipheno]*W_con + e_error
    
    sim_error[,ipheno] <- e_error
    sim_Y[,ipheno]     <- Y
  }
  
  # Test for differential co-expression using the scaled marginal framework
  SMAT_START = Sys.time()
  pval_SMAT_adj_dglm[isim] <- suppressWarnings(
    SMAT_cor_cov_dglm(sim_y = sim_Y, 
                      sim_g = X, 
                      sim_z = W_con,
                      includeXsquare = FALSE))
  SMAT_END = Sys.time()
  SMAT_temp <- 
    as.numeric(SMAT_END - SMAT_START, units="secs")
  SMAT_time = c(SMAT_time, SMAT_temp)
  
  # Test for differential co-expression using the scaled marginal framework
  #' Include the term X^2
  SMAT_includeXsquare_START = Sys.time()
  pval_SMAT_includeXsquare_adj_dglm[isim] <- suppressWarnings(
    SMAT_cor_cov_dglm(sim_y = sim_Y, 
                      sim_g = X, 
                      sim_z = W_con,
                      includeXsquare = TRUE))
  SMAT_includeXsquare_END = Sys.time()
  SMAT_includeXsquare_temp <- 
    as.numeric(SMAT_includeXsquare_END - SMAT_includeXsquare_START, units="secs")
  SMAT_includeXsquare_time = c(SMAT_includeXsquare_time, 
                               SMAT_includeXsquare_temp)
  
  
  # Test for differential co-expression using the Cauchy combination test
  CCT_START = Sys.time()
  pval_CCT_adj_dglm[isim] <- CCT_cor_cov_dglm(sim_y = sim_Y, 
                                              sim_g = X, 
                                              sim_z = W_con,
                                              includeXsquare = FALSE)
  CCT_END = Sys.time()
  CCT_temp <- 
    as.numeric(CCT_END - CCT_START, units="secs")
  CCT_time = c(CCT_time, CCT_temp)
  
  #' Test for differential co-expression using the Cauchy combination test
  #' Include the term X^2
  CCT_includeXsquare_START = Sys.time()
  pval_includeXsquare_CCT_adj_dglm[isim] <- CCT_cor_cov_dglm(sim_y = sim_Y, 
                                                             sim_g = X, 
                                                             sim_z = W_con,
                                                             includeXsquare = TRUE)
  CCT_includeXsquare_END = Sys.time()
  CCT_includeXsquare_temp <- 
    as.numeric(CCT_includeXsquare_END - CCT_includeXsquare_START, units="secs")
  CCT_includeXsquare_time = c(CCT_includeXsquare_time, 
                              CCT_includeXsquare_temp)
  
  end_time <- Sys.time()  
  one_sim_time <- as.numeric(end_time - start_time, units="secs")
  simulation_time = c(simulation_time, one_sim_time)
  
  if (isim %% 50 == 0) {
    print(isim)
  }
}

file_name_SMAT_pvalue               <- paste("pvalue", "SMAT", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
file_name_includeXsquar_SMAT_pvalue <- paste("pvalue", "includeXsquar_SMAT", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
file_name_CCT_pvalue                <- paste("pvalue", "CCT", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
file_name_includeXsquare_CCT_pvalue <- paste("pvalue", "includeXsquare_CCT", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")

file_SMAT_time                <- paste("time", "SMAT_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
file_includeXsquare_SMAT_time <- paste("time", "includeXsquare_SMAT_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
file_CCT_time                 <- paste("time", "CCT_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
file_includeXsquare_CCT_time  <- paste("time", "includeXsquare_CCT_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")
file_simulation_time          <- paste("time", "simulation_time", "adj", "dglm", Nsim, "JobSlice", JobSlice, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, sep = "_")

save(pval_SMAT_adj_dglm, 
     file = paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_name_SMAT_pvalue, ".RData", sep = ""))
save(pval_SMAT_includeXsquare_adj_dglm, 
     file = paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_name_includeXsquar_SMAT_pvalue, ".RData", sep = ""))
save(pval_CCT_adj_dglm, 
     file = paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_name_CCT_pvalue, ".RData", sep = ""))
save(pval_includeXsquare_CCT_adj_dglm, 
     file = paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_name_includeXsquare_CCT_pvalue, ".RData", sep = ""))

save(simulation_time, 
     file = paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_simulation_time, ".RData", sep = ""))
save(SMAT_time, 
     file = paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_SMAT_time, ".RData", sep = ""))
save(SMAT_includeXsquare_time, 
     file = paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_includeXsquare_SMAT_time, ".RData", sep = ""))
save(CCT_time, 
     file = paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_CCT_time, ".RData", sep = ""))
save(CCT_includeXsquare_time, 
     file = paste("Result/New_4_6_Simulation_ExcludeGAMUT/", file_includeXsquare_CCT_time, ".RData", sep = ""))

paste("----------End:", array.id, "out of", dim(param_comb)[1], "----------", sep = " ")
