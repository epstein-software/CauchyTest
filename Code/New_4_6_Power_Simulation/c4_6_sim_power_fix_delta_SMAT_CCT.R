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

Nsub_array   <- c(500, 1000, 10000, 30000, 50000)
Npheno_array <- c(4, 6, 8, 10, 12)
maf_array    <- 0.25
gamma_array  <- c(0.05, 0.1, 0.15, 0.2)
delta_array  <- seq(0.05, 0.75, by = 0.05)    
Nsim         <- 100000

# Nsub_array   <- c(500)
# Npheno_array <- c(4, 6)
# maf_array    <- c(0.15, 0.25)
# Nsim   <- 10

param_comb <- expand.grid(Nsub   = Nsub_array,
                          Npheno = Npheno_array,
                          maf    = maf_array,
                          gamma  = gamma_array,
                          delta  = delta_array)

Nsub   <- param_comb[array.id, "Nsub"]
Npheno <- param_comb[array.id, "Npheno"]
maf    <- param_comb[array.id, "maf"]
gamma  <- param_comb[array.id, "gamma"]
delta  <- param_comb[array.id, "delta"]

source('Code/New_4_6_Power_Simulation/c4_6_sim_source_code.R') 

#' Seed used for simulating the coefficients, the coefficients will be fixed in
#' the beginning. Instead of random simulating like this, we could 
int_array            <- array(NA,dim=c(Npheno))
main_beta_g_array    <- array(NA,dim=c(Npheno))
main_gamma_con_array <- array(NA,dim=c(Npheno))
main_delta_con_array <- array(NA,dim=c(Npheno))

for (sim_Npheno in 1:Npheno) {
  set.seed(20211101 + sim_Npheno)
  
  # randomly-generated intercept for phenotypes
  int_array[sim_Npheno] <- rnorm(1, 0, 5) 
  
  # randomly-generated main effect of SNP on phenotypes
  main_beta_g_array[sim_Npheno] <- runif(1, 0, 0.2)
  
  # randomly-generate main effect of continuous covariate on genotype
  main_gamma_con_array[sim_Npheno] <- gamma
  
  # randomly-generate main effect of interaction term
  main_delta_con_array[sim_Npheno] <- delta
}

#' Number of off-diagonal elements in covariance matrix of phenotype
#' (within subject)
Num_off_diag <- gamma(Npheno+1)/(gamma(3)*gamma(Npheno+1-2)) 

# Initialize arrays/vectors to contain p-values and time 
# pval_GAMuT_ProjMatrix_adj_dglm   <- array(NA,dim=c(Nsim))
# pval_GAMuT_LineKernel_adj_dglm   <- array(NA,dim=c(Nsim))
pval_SMAT_adj_dglm               <- array(NA,dim=c(Nsim))
pval_CCT_adj_dglm                <- array(NA,dim=c(Nsim))
pval_includeXsquare_CCT_adj_dglm <- array(NA,dim=c(Nsim))

# GAMuT_ProjMatrix_time = c()
# GAMuT_LineKernel_time = c()
SMAT_time = c()
CCT_time = c()
CCT_includeXsquare_time = c()
simulation_time = c()

for(isim in 1:Nsim) {
  
  start_time <- Sys.time()
  set.seed(2021 + isim)
  
  #' Each subject has their own X and W. Within each subject, X and W are fixed,
  #' but the coefficients vary again the phenotypes. The phenotypes vary by 
  #' randomly simulating the error term. This step results in (Nsub by Npheno)
  #'  simulated Y 
  X     <- rbinom(Nsub, 2, maf)
  W_con <- rnorm(Nsub, 0, 1)
  
  sim_Y     <- matrix(, nrow = Nsub, ncol = Npheno)
  sim_error <- matrix(, nrow = Nsub, ncol = Npheno)
  
  for (ipheno in 1:Npheno) {
    
    set.seed(2021 + isim + ipheno)
    
    e_error <- rnorm(Nsub, 0, 1)
    Y <- int_array[ipheno] + 
      main_beta_g_array[ipheno]*X + 
      main_gamma_con_array[ipheno]*W_con + 
      main_delta_con_array[ipheno]*X*W_con +
      e_error
    
    sim_error[,ipheno] <- e_error
    sim_Y[,ipheno]     <- Y
  }
  
  # Test for differential co-expression using the scaled marginal framework
  SMAT_START = Sys.time()
  pval_SMAT_adj_dglm[isim] <- suppressWarnings(
    SMAT_cor_cov_dglm(sim_y = sim_Y, 
                      sim_g = X, 
                      sim_z = W_con))
  SMAT_END = Sys.time()
  SMAT_time = c(SMAT_time, SMAT_END - SMAT_START)
  
  # Test for differential co-expression using the Cauchy combination test
  CCT_START = Sys.time()
  pval_CCT_adj_dglm[isim] <- CCT_cor_cov_dglm(sim_y = sim_Y, 
                                              sim_g = X, 
                                              sim_z = W_con,
                                              includeXsquare = FALSE)
  CCT_END = Sys.time()
  CCT_time = c(CCT_time, CCT_END - CCT_START)
  
  #' Test for differential co-expression using the Cauchy combination test
  #' Include the term X^2
  CCT_includeXsquare_START = Sys.time()
  pval_includeXsquare_CCT_adj_dglm[isim] <- CCT_cor_cov_dglm(sim_y = sim_Y, 
                                                             sim_g = X, 
                                                             sim_z = W_con,
                                                             includeXsquare = TRUE)
  CCT_includeXsquare_END = Sys.time()
  CCT_includeXsquare_time = c(CCT_includeXsquare_time, 
                              CCT_includeXsquare_END - CCT_includeXsquare_START)
  
  # # Test for differential co-expression using the GAMuT framework with LineKernel
  # GAMuT_LineKernel_START = Sys.time()
  # pval_GAMuT_LineKernel_adj_dglm[isim] <-
  #   GAMuT_cor_cov_dglm(sim_y = sim_Y,
  #                      sim_g = X,
  #                      sim_z = W_con,
  #                      kernel = "LineKernel")
  # GAMuT_LineKernel_END = Sys.time()
  # GAMuT_LineKernel_time = c(GAMuT_LineKernel_time, 
  #                           GAMuT_LineKernel_END - GAMuT_LineKernel_START)
  # 
  # GAMuT_ProjMatrix_START = Sys.time()
  # pval_GAMuT_ProjMatrix_adj_dglm[isim] <-
  #   GAMuT_cor_cov_dglm(sim_y = sim_Y,
  #                      sim_g = X,
  #                      sim_z = W_con,
  #                      kernel = "ProjMatrix")
  # GAMuT_ProjMatrix_END = Sys.time()
  # GAMuT_ProjMatrix_time = c(GAMuT_ProjMatrix_time, 
  #                           GAMuT_ProjMatrix_END - GAMuT_ProjMatrix_START)
  
  end_time <- Sys.time()  
  one_sim_time <- end_time - start_time
  simulation_time = c(simulation_time, one_sim_time)
  
  if (isim %% 100 == 0) {
    print(isim)
  }
}


# file_name_GAMuT_ProjMatrix_pvalue   <- paste("pvalue", "GAMuT_ProjMatrix", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, sep = "_")
# file_name_GAMuT_LineKernel_pvalue   <- paste("pvalue", "GAMuT_LineKernel", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, sep = "_")
file_name_SMAT_pvalue               <- paste("pvalue", "SMAT", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, sep = "_")
file_name_CCT_pvalue                <- paste("pvalue", "CCT", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, sep = "_")
file_name_includeXsquare_CCT_pvalue <- paste("pvalue", "includeXsquare_CCT", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, sep = "_")

# file_GAMuT_ProjMatrix_time   <- paste("time", "GAMuT_ProjMatrix_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, sep = "_")
# file_GAMuT_LineKernel_time   <- paste("time", "GAMuT_LineKernel_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, sep = "_")
file_SMAT_time               <- paste("time", "SMAT_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, sep = "_")
file_CCT_time                <- paste("time", "CCT_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, sep = "_")
file_includeXsquare_CCT_time <- paste("time", "includeXsquare_CCT_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, sep = "_")
file_simulation_time         <- paste("time", "simulation_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, sep = "_")


# save(pval_GAMuT_ProjMatrix_adj_dglm, 
#      file = paste("Result/New_4_6_Power_Simulation/", file_name_GAMuT_ProjMatrix_pvalue, ".RData", sep = ""))
# save(pval_GAMuT_LineKernel_adj_dglm, 
#      file = paste("Result/New_4_6_Power_Simulation/", file_name_GAMuT_LineKernel_pvalue, ".RData", sep = ""))
save(pval_SMAT_adj_dglm, 
     file = paste("Result/New_4_6_Power_Simulation/", file_name_SMAT_pvalue, ".RData", sep = ""))
save(pval_CCT_adj_dglm, 
     file = paste("Result/New_4_6_Power_Simulation/", file_name_CCT_pvalue, ".RData", sep = ""))
save(pval_includeXsquare_CCT_adj_dglm, 
     file = paste("Result/New_4_6_Power_Simulation/", file_name_includeXsquare_CCT_pvalue, ".RData", sep = ""))

save(simulation_time, 
     file = paste("Result/New_4_6_Power_Simulation/", file_simulation_time, ".RData", sep = ""))
# save(GAMuT_ProjMatrix_time, 
#      file = paste("Result/New_4_6_Power_Simulation/", file_GAMuT_ProjMatrix_time, ".RData", sep = ""))
# save(GAMuT_LineKernel_time, 
#      file = paste("Result/New_4_6_Power_Simulation/", file_GAMuT_LineKernel_time, ".RData", sep = ""))
save(SMAT_time, 
     file = paste("Result/New_4_6_Power_Simulation/", file_SMAT_time, ".RData", sep = ""))
save(CCT_time, 
     file = paste("Result/New_4_6_Power_Simulation/", file_CCT_time, ".RData", sep = ""))
save(CCT_includeXsquare_time, 
     file = paste("Result/New_4_6_Power_Simulation/", file_includeXsquare_CCT_time, ".RData", sep = ""))


