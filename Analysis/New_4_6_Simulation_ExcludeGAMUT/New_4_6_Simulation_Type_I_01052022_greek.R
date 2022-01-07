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

library(tidyverse)
# library(qqman)
library(lattice)

# setwd("/n/holystore01/LABS/xlin/Lab/xihaoli/DOHG/")
setwd("~/Dropbox/Emory Courses/DOHG/")
source("Analysis/New_4_6_Simulation_ExcludeGAMUT/qqplot-source-code.R")

###### Predefined global parameter
#' Number of simulations, number of subjects, Number of phenotypes and MAF
#' (independent of covariates)
Nsub_array   <- c(5000, 10000)
Npheno_array <- c(4, 6, 8, 10, 12)
maf_array    <- c(0.25)
gamma_array  <- c(NA, 0.15, 0.55, 0.95)

param_comb <- expand.grid(Nsub   = Nsub_array,
                          gamma  = gamma_array,
                          maf    = maf_array,
                          Npheno = Npheno_array)

###### Plot 5,000 and 10,000 subjects on the same plot given gamma and Npheno
rowN <- nrow(param_comb)/2

for (index in 1:rowN) {
  
  # Nsub = 5,000
  slice_inde_1 = param_comb[index*2-1,]
  file_combine_pvalue_1 <- paste("pvalue", "combine", "adj", "dglm", "subj", slice_inde_1["Nsub"], "pheno", slice_inde_1["Npheno"], "MAF", slice_inde_1["maf"], "Gamma", slice_inde_1["gamma"], sep = "_")
  assign(file_combine_pvalue_1, get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT_Summary/pval_combine/", file_combine_pvalue_1, ".RData", sep = ""))))
  dat_1 <- get(file_combine_pvalue_1)
  
  my.pvalue.list.1 <-list("CCT"                 = dat_1$pval_CCT_adj_dglm_combine,
                          "CCT_includeXsquare"  = dat_1$pval_includeXsquare_CCT_adj_dglm_combine,
                          "SMAT"                = dat_1$pval_SMAT_adj_dglm_combine, 
                          "SMAT_includeXsquare" = dat_1$pval_SMAT_includeXsquare_adj_dglm_combine)
  
  lambda1_CCT                 <- dchisq(qchisq(median(dat_1$pval_CCT_adj_dglm_combine),df=1),df=1)/dchisq(qchisq(0.5,df=1),df=1)
  lambda1_CCT_includeXsquare  <- dchisq(qchisq(median(dat_1$pval_includeXsquare_CCT_adj_dglm_combine),df=1),df=1)/dchisq(qchisq(0.5,df=1),df=1)
  lambda1_SMAT                <- dchisq(qchisq(median(dat_1$pval_SMAT_adj_dglm_combine),df=1),df=1)/dchisq(qchisq(0.5,df=1),df=1)
  lambda1_SMAT_includeXsquare <- dchisq(qchisq(median(dat_1$pval_SMAT_includeXsquare_adj_dglm_combine),df=1),df=1)/dchisq(qchisq(0.5,df=1),df=1)
  
  p1 <- qqunif.plot(my.pvalue.list.1, 
              auto.key=list(corner=c(.95,.05)),
              col = c("#003B73", "#60A3D9", "#778A35", "#ECF87F"),
              key=list(space="right",
                       lines=list(col= c("#003B73", "#60A3D9", "#778A35", "#ECF87F"), 
                                  lwd=6),
                       text=list(
                         c(paste("CCT,                       \u03BB = ", round(lambda1_CCT, 3)),
                           paste("CCT include X^2,   \u03BB = ", round(lambda1_CCT_includeXsquare, 3)),
                           paste("SMAT,                    \u03BB = ", round(lambda1_SMAT, 3)),
                           paste("SMAT include X^2, \u03BB = ", round(lambda1_SMAT_includeXsquare, 3))))),
              main="5,000 Subjects")
  
  # Nsub = 10,000
  slice_inde_2 = param_comb[index*2,]
  file_combine_pvalue_2 <- paste("pvalue", "combine", "adj", "dglm", "subj", slice_inde_2["Nsub"], "pheno", slice_inde_2["Npheno"], "MAF", slice_inde_2["maf"], "Gamma", slice_inde_2["gamma"], sep = "_")
  assign(file_combine_pvalue_2, get(load(paste("Result/New_4_6_Simulation_ExcludeGAMUT_Summary/pval_combine/", file_combine_pvalue_2, ".RData", sep = ""))))
  dat_2 <- get(file_combine_pvalue_2)
  
  my.pvalue.list.2 <-list("CCT"                 = dat_2$pval_CCT_adj_dglm_combine,
                          "CCT_includeXsquare"  = dat_2$pval_includeXsquare_CCT_adj_dglm_combine,
                          "SMAT"                = dat_2$pval_SMAT_adj_dglm_combine, 
                          "SMAT_includeXsquare" = dat_2$pval_SMAT_includeXsquare_adj_dglm_combine)
  
  lambda2_CCT                 <- dchisq(qchisq(median(dat_2$pval_CCT_adj_dglm_combine),df=1),df=1)/dchisq(qchisq(0.5,df=1),df=1)
  lambda2_CCT_includeXsquare  <- dchisq(qchisq(median(dat_2$pval_includeXsquare_CCT_adj_dglm_combine),df=1),df=1)/dchisq(qchisq(0.5,df=1),df=1)
  lambda2_SMAT                <- dchisq(qchisq(median(dat_2$pval_SMAT_adj_dglm_combine),df=1),df=1)/dchisq(qchisq(0.5,df=1),df=1)
  lambda2_SMAT_includeXsquare <- dchisq(qchisq(median(dat_2$pval_SMAT_includeXsquare_adj_dglm_combine),df=1),df=1)/dchisq(qchisq(0.5,df=1),df=1)
  
  p2 <- qqunif.plot(my.pvalue.list.2, 
              auto.key=list(corner=c(.95,.05)),
              col = c("#003B73", "#60A3D9", "#778A35", "#ECF87F"),
              key=list(space="right",
                       lines=list(col= c("#003B73", "#60A3D9", "#778A35", "#ECF87F"), 
                                  lwd=6),
                       text=list(
                         c(paste("CCT,                      \u03BB = ", round(lambda2_CCT, 3)),
                           paste("CCT include X^2,   \u03BB = ", round(lambda2_CCT_includeXsquare, 3)),
                           paste("SMAT,                    \u03BB = ", round(lambda2_SMAT, 3)),
                           paste("SMAT include X^2, \u03BB = ", round(lambda2_SMAT_includeXsquare, 3))))),
              main="10,000 Subjects")
  
  if (is.na(slice_inde_2["gamma"][[1]])){
    gamma_text <- "uniform(0, 0.2)"
  } else {
    gamma_text <- slice_inde_2["gamma"]
  }
  
  file_name = paste("Analysis/New_4_6_Simulation_ExcludeGAMUT/type_i_plot_01052022/",
                    "type_I_error_index_", index, "_pheno_", slice_inde_1["Npheno"], "_Gamma_", slice_inde_1["gamma"], ".png", sep = "")
  png(file_name, width = 4800,height=6000, res = 600)
  grid.arrange(p1,p2, ncol=1, nrow=2,
               top = textGrob(paste("Type I error p-value: ", "# of phenotype is ", slice_inde_1["Npheno"], ", \u03B3 =", gamma_text, sep = ""),
                              gp=gpar(fontsize=14,font=3)))
  dev.off()
}
