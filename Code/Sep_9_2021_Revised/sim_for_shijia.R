############################### Install packages ###############################
# install.packages("MASS")
# install.packages("gdata")
# install.packages("CompQuadForm")
# install.packages("dglm")

library(MASS)
library(gdata)
library(CompQuadForm)
library(dglm)

setwd("C:/Users/sbian6/Dropbox/Emory Courses/DOHG")

# source code for the DKAT method
# source('./Code/Sep_9_2021_Revised/DKAT_source.R')
# source code for scaled marginal model
source('./Code/Sep_9_2021_Revised/SMAT_code.R') 
# source code for remaining functions
# source('Code/Sep_9_2021_Revised/source_files_for_taylor.R') 
source('Code/Sep_9_2021_Revised/source_files_for_shijia.R')

################################################################################
set.seed(500)

Nsim <- 10000 # Number of simulations

Nsub <- 500 # Number of subjects
Ngenes <- 4 # Number of genes in pathway
maf <- 0.3 # Minor allele frequency of test genotype that doesn't depend on covariate
maf_0 <- 0.3 # minor-allele frequency of test genotype when binary covariate has value of 0
maf_1 <- 0.8 # minor-allele frequency of test genotype when binary covariate has value of 1
p_z_bin <- 0.5
beta_z_bin <- 0
beta_z_con <- 0.2

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

#Initialize arrays/vectors to contain p-values

pval_GAMuT_adj_dglm <- array(NA,dim=c(Nsim,2))
pval_SMAT_adj_dglm <- array(NA,dim=c(Nsim))
pval_CCT_adj_dglm <- array(NA,dim=c(Nsim))

for(isim in 1:Nsim) {
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
  
  # Test for differential co-expression using the GAMuT framework
  pval_GAMuT_adj_dglm[isim, ] <-
    GAMuT_cor_cov_dglm(sim_y, sim_g, sim_z_con)
  
  # Test for differential co-expression using the scaled marginal framework
  pval_SMAT_adj_dglm[isim] <- SMAT_cor_cov_dglm(sim_y, sim_g, sim_z_con)
  
  # Test for differential co-expression using the Cauchy combination test
  pval_CCT_adj_dglm[isim] <- CCT_cor_cov_dglm(sim_y, sim_g, sim_z_con)
  
  if (isim %% 100 == 0) {
    print(isim)
  }
  
}

hist(pval_SMAT_adj_dglm)
hist(pval_CCT_adj_dglm)

save(pval_GAMuT_adj_dglm, file = "Result/Sep132021/pval_GAMuT_adj_dglm_shijia.RData")
save(pval_SMAT_adj_dglm, file = "Result/Sep132021/pval_SMAT_adj_dglm_shijia.RData")
save(pval_CCT_adj_dglm, file = "Result/Sep132021/pval_CCT_adj_dglm_shijia.RData")


load("Result/Sep132021/pval_GAMuT_adj_dglm_shijia.RData")
load("Result/Sep132021/pval_SMAT_adj_dglm_shijia.RData")
load("Result/Sep132021/pval_CCT_adj_dglm_shijia.RData")


hist(pval_SMAT_adj_dglm)
hist(pval_CCT_adj_dglm)

sum(pval_SMAT_adj_dglm < 0.05)/length(pval_SMAT_adj_dglm)
# [1] 0.0535
sum(pval_CCT_adj_dglm < 0.05)/length(pval_CCT_adj_dglm)
# [1] 0.0481

