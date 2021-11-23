############################### Install packages ###############################
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
# setwd("C:/Users/sbian6/Dropbox/Emory Courses/DOHG")
# source code for scaled marginal model
source('Code/c4_6_Simulation_HPC/SMAT_Code.R') 
# source code for remaining functions
source('Code/c4_6_Simulation_HPC/source_files_for_shijia.R')

start_time <- Sys.time()
################################################################################

Nsim <- 10000 # Number of simulations
Nsub <- 1000 # Number of subjects
Ngenes <- 16 # Number of genes in pathway
maf <- 0.3 # Minor allele frequency of test genotype that doesn't depend on covariate
maf_0 <- 0.3 # minor-allele frequency of test genotype when binary covariate has value of 0
maf_1 <- 0.8 # minor-allele frequency of test genotype when binary covariate has value of 1
p_z_bin <- 0.5
beta_z_bin <- 0
beta_z_con <- 0.2

#Initialize arrays/vectors to contain p-values

pval_GAMuT_adj_dglm <- array(NA,dim=c(Nsim,2))
pval_SMAT_adj_dglm <- array(NA,dim=c(Nsim))
pval_CCT_adj_dglm <- array(NA,dim=c(Nsim))

GAMuT_time = c()
SMAT_time = c()
CCT_time = c()
for(isim in 1:Nsim) {
  set.seed(2021123 + isim)
  
  int <- rnorm(Ngenes, 0, 5) # randomly-generated intercept for gene expression
  main_beta_g <- runif(Ngenes, 0, 0.2) # randomly-generated main effect of SNP on gene expression
  var_beta_g <- runif(Ngenes, 0, 0.2) # random-generated variance effect of SNP on gene expression
  
  main_beta_z_bin <- runif(Ngenes,0,beta_z_bin) # randomly-generated main effect of binary covariate on gene expression
  main_beta_z_con <- runif(Ngenes,0,beta_z_con) # randomly-generate main effect of continuous covariate on gene expression
  var_beta_z_bin <- runif(Ngenes,0,beta_z_bin) # randomly-generated variance effect of binary covariate on gene expression
  var_beta_z_con <- runif(Ngenes,0,beta_z_con) # randomly-generated variance effect of continuous covariate on gene expression
  
  Num_off_diag <- gamma(Ngenes+1)/(gamma(3)*gamma(Ngenes+1-2)) # Number of off-diagonal elements in covariance matrix of gene expression (within subject)
  
  cov_beta_g <- runif(Num_off_diag,0,0) # randomly generated off-diagnonal covariance of gene expression (as function of SNP)
  cov_beta_z_bin <- runif(Num_off_diag,0,0)
  cov_beta_z_con <- runif(Num_off_diag,0,0)
  
  sigma_e <- runif(Ngenes, 1, 1) # randomly generated non-genetic variance of gene expression
  
  # generates binary environmental covariate
  sim_z_bin <- rbinom(Nsub, 1, p_z_bin) 
  
  sim_z_con <- rnorm(Nsub, 0, 1)
  
  # generate SNP genotype that is not dependent on binary covariate
  sim_g <- rbinom(Nsub, 2, maf) 
  
  preds2 <- rbind(sim_g, 
                  sim_z_bin, 
                  sim_z_con)
  preds3 <- split(preds2, 
                  rep(1:ncol(preds2), each = nrow(preds2)))
  
  #Generate expression data
  sim_y <-
    t(
      sapply(
        preds3,
        gen_express_cov,
        int = int,
        main_beta_g = main_beta_g,
        var_beta_g = var_beta_g,
        Num_off_diag = Num_off_diag,
        cov_beta_g = cov_beta_g,
        sigma_e = sigma_e,
        main_beta_z_bin = main_beta_z_bin,
        main_beta_z_con = main_beta_z_con,
        var_beta_z_bin = var_beta_z_bin,
        var_beta_z_con = var_beta_z_con,
        cov_beta_z_bin = cov_beta_z_bin,
        cov_beta_z_con = cov_beta_z_con
      )
    )
  
  GAMuT_START = Sys.time()
  # Test for differential co-expression using the GAMuT framework
  pval_GAMuT_adj_dglm[isim, ] <-
    GAMuT_cor_cov_dglm(sim_y, sim_g, sim_z_con)
  GAMuT_END = Sys.time()
  GAMuT_time = c(GAMuT_time, GAMuT_END - GAMuT_START)
  
  SMAT_START = Sys.time()
  # Test for differential co-expression using the scaled marginal framework
  pval_SMAT_adj_dglm[isim] <- SMAT_cor_cov_dglm(sim_y, sim_g, sim_z_con)
  SMAT_END = Sys.time()
  SMAT_time = c(SMAT_time, SMAT_END - SMAT_START)
  
  CCT_START = Sys.time()
  # Test for differential co-expression using the Cauchy combination test
  pval_CCT_adj_dglm[isim] <- CCT_cor_cov_dglm(sim_y, sim_g, sim_z_con)
  CCT_END = Sys.time()
  CCT_time = c(CCT_time, CCT_END - CCT_START)
  
  if (isim %% 100 == 0) {
    print(isim)
  }
  
}

end_time <- Sys.time()

file_name_GAMuT <- paste(Sys.Date(), "pvalue", "GAMuT", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Ngenes, sep = "_")
file_name_SMAT <- paste(Sys.Date(), "pvalue", "SMAT", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Ngenes, sep = "_")
file_name_CCT <- paste(Sys.Date(), "pvalue", "CCT", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Ngenes, sep = "_")
file_time_begin <- paste(Sys.Date(), "pvalue", "time_begin", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Ngenes, sep = "_")
file_time_end <- paste(Sys.Date(), "pvalue", "time_end", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Ngenes, sep = "_")
file_GAMuT_time <- paste(Sys.Date(), "pvalue", "GAMuT_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Ngenes, sep = "_")
file_CCT_time <- paste(Sys.Date(), "pvalue", "CCT_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Ngenes, sep = "_")
file_SMAT_time <- paste(Sys.Date(), "pvalue", "SMAT_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Ngenes, sep = "_")


save(pval_GAMuT_adj_dglm, 
     file = paste("Result/c4_6_Simulation_HPC/", file_name_GAMuT, ".RData", sep = ""))
save(pval_SMAT_adj_dglm, 
     file = paste("Result/c4_6_Simulation_HPC/", file_name_SMAT, ".RData", sep = ""))
save(pval_CCT_adj_dglm, 
     file = paste("Result/c4_6_Simulation_HPC/", file_name_CCT, ".RData", sep = ""))
save(start_time, 
     file = paste("Result/c4_6_Simulation_HPC/", file_time_begin, ".RData", sep = ""))
save(end_time, 
     file = paste("Result/c4_6_Simulation_HPC/", file_time_end, ".RData", sep = ""))
save(GAMuT_time, 
     file = paste("Result/c4_6_Simulation_HPC/", file_GAMuT_time, ".RData", sep = ""))
save(CCT_time, 
     file = paste("Result/c4_6_Simulation_HPC/", file_CCT_time, ".RData", sep = ""))
save(SMAT_time, 
     file = paste("Result/c4_6_Simulation_HPC/", file_SMAT_time, ".RData", sep = ""))

# print(end_time - start_time)

# load("Result/Sep132021/pval_GAMuT_adj_dglm_shijia.RData")
# load("Result/Sep132021/pval_SMAT_adj_dglm_shijia.RData")
# load("Result/Sep132021/pval_CCT_adj_dglm_shijia.RData")
# 
# 
# hist(pval_SMAT_adj_dglm)
# hist(pval_CCT_adj_dglm)
# 
# sum(pval_SMAT_adj_dglm < 0.05)/length(pval_SMAT_adj_dglm)
# # [1] 0.0535
# sum(pval_CCT_adj_dglm < 0.05)/length(pval_CCT_adj_dglm)
# # [1] 0.0481
# 
