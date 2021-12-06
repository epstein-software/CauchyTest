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

setwd("~/Dropbox/Emory Courses/DOHG/")

###### Predefined global parameter
#' Number of simulations, number of subjects, Number of phenotypes and MAF
#' (independent of covariates)

Nsub_array   <- c(500)
Npheno_array <- c(4, 6, 8, 10, 12)
maf_array    <- 0.25
gamma_array  <- c(0.15, 0.25, 0.35)
delta_array  <- c(0.5, 0.75, 1) 
delta_select_array <- c(0.25, 0.5, 0.75, 1)
Nsim         <- 1000

param_comb <- expand.grid(Nsub   = Nsub_array,
                          Npheno = Npheno_array,
                          maf    = maf_array,
                          gamma  = gamma_array,
                          delta  = delta_array,
                          delta_select = delta_select_array)

############ ---------- Construct the Overall Dataframe ----------- ############

full_df <- data.frame()
for (array.id in 1:nrow(param_comb)) {
  
  Nsub   <- param_comb[array.id, "Nsub"]
  Npheno <- param_comb[array.id, "Npheno"]
  maf    <- param_comb[array.id, "maf"]
  gamma  <- param_comb[array.id, "gamma"]
  delta  <- param_comb[array.id, "delta"]
  delta_select  <- param_comb[array.id, "delta_select"]
  
  file_SMAT_time               <- paste("time", "SMAT_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_CCT_time                <- paste("time", "CCT_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_includeXsquare_CCT_time <- paste("time", "includeXsquare_CCT_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_simulation_time         <- paste("time", "simulation_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  
  file_name_SMAT_pvalue               <- paste("pvalue", "SMAT", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_name_CCT_pvalue                <- paste("pvalue", "CCT", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_name_includeXsquare_CCT_pvalue <- paste("pvalue", "includeXsquare_CCT", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")

  assign(file_SMAT_time,
         get(load(paste("Result/New_4_6_Power_Simulation/", file_SMAT_time, ".RData", sep = ""))))
  assign(file_CCT_time,
         get(load(paste("Result/New_4_6_Power_Simulation/", file_CCT_time, ".RData", sep = ""))))
  assign(file_includeXsquare_CCT_time,
         get(load(paste("Result/New_4_6_Power_Simulation/", file_includeXsquare_CCT_time, ".RData", sep = ""))))
  assign(file_simulation_time,
         get(load(paste("Result/New_4_6_Power_Simulation/", file_simulation_time, ".RData", sep = ""))))

  assign(file_name_SMAT_pvalue,
         get(load(paste("Result/New_4_6_Power_Simulation/", file_name_SMAT_pvalue, ".RData", sep = ""))))
  assign(file_name_CCT_pvalue,
         get(load(paste("Result/New_4_6_Power_Simulation/", file_name_CCT_pvalue, ".RData", sep = ""))))
  assign(file_name_includeXsquare_CCT_pvalue,
         get(load(paste("Result/New_4_6_Power_Simulation/", file_name_includeXsquare_CCT_pvalue, ".RData", sep = ""))))

  temp_DF <- data.frame(
    "Nsub"          = Nsub,
    "Npheno"        = Npheno,
    "maf"           = maf,
    "gamma"         = gamma,
    "delta"         = delta,
    "delta_select"  = delta_select,
    "SMAT"             = get(load(paste("Result/New_4_6_Power_Simulation/", file_name_SMAT_pvalue, ".RData", sep = ""))),
    "CCT"              = get(load(paste("Result/New_4_6_Power_Simulation/", file_name_CCT_pvalue, ".RData", sep = ""))),
    "CCT_includeXsquare" = get(load(paste("Result/New_4_6_Power_Simulation/", file_name_includeXsquare_CCT_pvalue, ".RData", sep = "")))
#    "SMAT_Time"             = get(load(paste("Result/New_4_6_Power_Simulation/", file_SMAT_time, ".RData", sep = ""))),
#    "CCT_Time"              = get(load(paste("Result/New_4_6_Power_Simulation/", file_CCT_time, ".RData", sep = ""))),
#    "CCT_includeXsquare_Time" = get(load(paste("Result/New_4_6_Power_Simulation/", file_includeXsquare_CCT_time, ".RData", sep = ""))),
#    "Overall_Time" = get(load(paste("Result/New_4_6_Power_Simulation/", file_simulation_time, ".RData", sep = "")))
  )
  
  full_df <- rbind(full_df, temp_DF)
}

#### ---- At Different Combinations of Number of Phenotypes and Delta ---- #####

#### ----                       Phenotype == 12                       ---- #####

######### Delta == 1
##--- Phenotype == 12, Delta == 1, delta_select_array == 1
dat <- full_df %>%
  filter(Npheno == 12,
         delta == 1,
         delta_select == 1) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")
  

##--- Phenotype == 12, Delta == 1, delta_select_array == 0.75
dat <- full_df %>%
  filter(Npheno == 12,
         delta == 1,
         delta_select == 0.75) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

##--- Phenotype == 12, Delta == 1, delta_select_array == 0.5
dat <- full_df %>%
  filter(Npheno == 12,
         delta == 1,
         delta_select == 0.5) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


##--- Phenotype == 12, Delta == 1, delta_select_array == 0.25
dat <- full_df %>%
  filter(Npheno == 12,
         delta == 1,
         delta_select == 0.25) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

######### Delta == 0.5
##--- Phenotype == 12, Delta == 0.5, delta_select_array == 1
dat <- full_df %>%
  filter(Npheno == 12,
         delta == 0.5,
         delta_select == 1) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


##--- Phenotype == 12, Delta == 0.5, delta_select_array == 0.75
dat <- full_df %>%
  filter(Npheno == 12,
         delta == 0.5,
         delta_select == 0.75) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

##--- Phenotype == 12, Delta == 0.5, delta_select_array == 0.5
dat <- full_df %>%
  filter(Npheno == 12,
         delta == 0.5,
         delta_select == 0.5) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


##--- Phenotype == 12, Delta == 0.5, delta_select_array == 0.25
dat <- full_df %>%
  filter(Npheno == 12,
         delta == 0.5,
         delta_select == 0.25) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


##--- Phenotype == 8, Delta == 0.75
dat <- full_df %>%
  filter(Npheno == 8,
         delta == 0.75) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

##--- Phenotype == 4, Delta == 0.75
dat <- full_df %>%
  filter(Npheno == 4,
         delta == 0.75) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

##--- Phenotype == 12, Delta == 0.5
dat <- full_df %>%
  filter(Npheno == 12,
         delta == 0.5) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

##--- Phenotype == 8, Delta == 0.5
dat <- full_df %>%
  filter(Npheno == 8,
         delta == 0.5) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


##--- Phenotype == 4, Delta == 0.5
dat <- full_df %>%
  filter(Npheno == 4,
         delta == 0.5) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

#### ----                       Phenotype == 8                       ---- #####

######### Delta == 1
##--- Phenotype == 8, Delta == 1, delta_select_array == 1
dat <- full_df %>%
  filter(Npheno == 8,
         delta == 1,
         delta_select == 1) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


##--- Phenotype == 8, Delta == 1, delta_select_array == 0.75
dat <- full_df %>%
  filter(Npheno == 8,
         delta == 1,
         delta_select == 0.75) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

##--- Phenotype == 8, Delta == 1, delta_select_array == 0.5
dat <- full_df %>%
  filter(Npheno == 8,
         delta == 1,
         delta_select == 0.5) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


##--- Phenotype == 8, Delta == 1, delta_select_array == 0.25
dat <- full_df %>%
  filter(Npheno == 8,
         delta == 1,
         delta_select == 0.25) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

######### Delta == 0.5
##--- Phenotype == 8, Delta == 0.5, delta_select_array == 1
dat <- full_df %>%
  filter(Npheno == 8,
         delta == 0.5,
         delta_select == 1) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


##--- Phenotype == 8, Delta == 0.5, delta_select_array == 0.75
dat <- full_df %>%
  filter(Npheno == 8,
         delta == 0.5,
         delta_select == 0.75) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

##--- Phenotype == 8, Delta == 0.5, delta_select_array == 0.5
dat <- full_df %>%
  filter(Npheno == 8,
         delta == 0.5,
         delta_select == 0.5) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


##--- Phenotype == 8, Delta == 0.5, delta_select_array == 0.25
dat <- full_df %>%
  filter(Npheno == 8,
         delta == 0.5,
         delta_select == 0.25) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


#### ----                       Phenotype == 4                       ---- #####

######### Delta == 1
##--- Phenotype == 4, Delta == 1, delta_select_array == 1
dat <- full_df %>%
  filter(Npheno == 4,
         delta == 1,
         delta_select == 1) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


##--- Phenotype == 4, Delta == 1, delta_select_array == 0.75
dat <- full_df %>%
  filter(Npheno == 4,
         delta == 1,
         delta_select == 0.75) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

##--- Phenotype == 4, Delta == 1, delta_select_array == 0.5
dat <- full_df %>%
  filter(Npheno == 4,
         delta == 1,
         delta_select == 0.5) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


##--- Phenotype == 4, Delta == 1, delta_select_array == 0.25
dat <- full_df %>%
  filter(Npheno == 4,
         delta == 1,
         delta_select == 0.25) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

######### Delta == 0.5
##--- Phenotype == 4, Delta == 0.5, delta_select_array == 1
dat <- full_df %>%
  filter(Npheno == 4,
         delta == 0.5,
         delta_select == 1) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


##--- Phenotype == 4, Delta == 0.5, delta_select_array == 0.75
dat <- full_df %>%
  filter(Npheno == 4,
         delta == 0.5,
         delta_select == 0.75) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

##--- Phenotype == 4, Delta == 0.5, delta_select_array == 0.5
dat <- full_df %>%
  filter(Npheno == 4,
         delta == 0.5,
         delta_select == 0.5) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")


##--- Phenotype == 4, Delta == 0.5, delta_select_array == 0.25
dat <- full_df %>%
  filter(Npheno == 4,
         delta == 0.5,
         delta_select == 0.25) %>%
  select(gamma, SMAT, CCT, CCT_includeXsquare) %>%
  gather(., "method", "pvalue", SMAT:CCT_includeXsquare) %>%
  group_by(gamma, method) %>%
  summarise(count_less_005       = sum(pvalue <= 0.05), 
            count_less_005_10000 = sum(pvalue <= 0.05/10000),
            count_less_005_20000 = sum(pvalue <= 0.05/20000),
            count_less_005_30000 = sum(pvalue <= 0.05/30000),
            count_less_005_40000 = sum(pvalue <= 0.05/40000),
            count_less_005_50000 = sum(pvalue <= 0.05/50000),
            count_less_005_60000 = sum(pvalue <= 0.05/60000),
            .groups="drop")

#' Given the same number of phenotype, maf (same for all originally), 
#' delta, delta_select and method, p-value are different for different gamma
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

######### Phenotype is 8
gamma_015 <- 
  full_df %>%
  filter(Npheno == 8,
         delta == 1,
         delta_select == 0.75,
         gamma == 0.15)
sum(gamma_015$SMAT <= 0.05/50000) # 7
specify_decimal(mean(gamma_015$SMAT), 20) # 0.01343130059667276693
specify_decimal(mean(gamma_015$CCT), 20)  # 0.00043311829323283403
specify_decimal(mean(gamma_015$CCT_includeXsquare), 20)  # 0.00034450325588403154

specify_decimal(median(gamma_015$SMAT), 20) # 0.00247922534895639313
specify_decimal(median(gamma_015$CCT), 20)  # 0.00001122445
specify_decimal(median(gamma_015$CCT_includeXsquare), 20)  # 0.0000008520913

specify_decimal(max(gamma_015$SMAT), 20) # 0.48949447667456835731
specify_decimal(max(gamma_015$CCT), 20)  # 0.02737312842337480401
specify_decimal(max(gamma_015$CCT_includeXsquare), 20)  # 0.07169513212008116199

specify_decimal(min(gamma_015$SMAT), 20) # 5.690977e-09
specify_decimal(min(gamma_015$CCT), 20)  # 3.183187e-11
specify_decimal(min(gamma_015$CCT_includeXsquare), 20)  # 4.4409e-16

gamma_025 <- 
  full_df %>%
  filter(Npheno == 8,
         delta == 1,
         delta_select == 0.75,
         gamma == 0.25)
sum(gamma_025$SMAT <= 0.05/50000) # 7
specify_decimal(mean(gamma_025$SMAT), 20) # 0.01343130059667275132
specify_decimal(mean(gamma_025$CCT), 20)  # 0.00043311829323283468
specify_decimal(mean(gamma_025$CCT_includeXsquare), 20)  # 0.00034450325588407350

specify_decimal(median(gamma_025$SMAT), 20) # 0.00247922534895639313
specify_decimal(median(gamma_025$CCT), 20)  # 0.00001122445
specify_decimal(median(gamma_025$CCT_includeXsquare), 20)  # 0.0000008520913

specify_decimal(max(gamma_025$SMAT), 20) # 0.48949447667455991962
specify_decimal(max(gamma_025$CCT), 20)  # 0.02737312842337469299
specify_decimal(max(gamma_025$CCT_includeXsquare), 20)  # 0.07169513212008526981

specify_decimal(min(gamma_025$SMAT), 20) # 5.690977e-09
specify_decimal(min(gamma_025$CCT), 20)  # 3.183187e-11
specify_decimal(min(gamma_025$CCT_includeXsquare), 20)  # 4.4409e-16

gamma_035 <- 
  full_df %>%
  filter(Npheno == 8,
         delta == 1,
         delta_select == 0.75,
         gamma == 0.35)
sum(gamma_035$SMAT <= 0.05/50000) # 7
specify_decimal(mean(gamma_035$SMAT), 20)  # 0.01343130059667275999
specify_decimal(mean(gamma_035$CCT), 20)   # 0.00043311829323283712
specify_decimal(mean(gamma_035$CCT_includeXsquare), 20)  # 0.00034450325588423938

specify_decimal(median(gamma_035$SMAT), 20) # 0.00247922534895639313
specify_decimal(median(gamma_035$CCT), 20)  # 0.00001122445
specify_decimal(median(gamma_035$CCT_includeXsquare), 20)  # 0.0000008520913

specify_decimal(max(gamma_035$SMAT), 20) # 0.48949447667456946753
specify_decimal(max(gamma_035$CCT), 20)  # 0.02737312842337469299
specify_decimal(max(gamma_035$CCT_includeXsquare), 20)  # 0.07169513212008515879

specify_decimal(min(gamma_035$SMAT), 20) # 5.690977e-09
specify_decimal(min(gamma_035$CCT), 20)  # 3.183187e-11
specify_decimal(min(gamma_035$CCT_includeXsquare), 20)  # 4.4409e-16

######### Phenotype is 12
gamma_015 <- 
  full_df %>%
  filter(Npheno == 12,
         delta == 1,
         delta_select == 0.75,
         gamma == 0.15)
sum(gamma_015$SMAT <= 0.05/50000) # 13
specify_decimal(mean(gamma_015$SMAT), 20) # 0.02245105569822626254
specify_decimal(mean(gamma_015$CCT), 20)  # 0.00017748158796399526
specify_decimal(mean(gamma_015$CCT_includeXsquare), 20)  # 7.927125e-05

specify_decimal(median(gamma_015$SMAT), 20) # 0.00275967613626793007
specify_decimal(median(gamma_015$CCT), 20)  # 0.00001122445
specify_decimal(median(gamma_015$CCT_includeXsquare), 20)  # 0.0000008520913

specify_decimal(max(gamma_015$SMAT), 20) # 0.78691592406092625289
specify_decimal(max(gamma_015$CCT), 20)  # 0.02158895824613304981
specify_decimal(max(gamma_015$CCT_includeXsquare), 20)  # 0.00913191681141212186

specify_decimal(min(gamma_015$SMAT), 20) # 8.618284e-11
specify_decimal(min(gamma_015$CCT), 20)  # 1.044564e-11
specify_decimal(min(gamma_015$CCT_includeXsquare), 20)  # 8.8818e-16

gamma_025 <- 
  full_df %>%
  filter(Npheno == 12,
         delta == 1,
         delta_select == 0.75,
         gamma == 0.25)
sum(gamma_025$SMAT <= 0.05/50000) # 13
specify_decimal(mean(gamma_025$SMAT), 20) # 0.02245105569822614111
specify_decimal(mean(gamma_025$CCT), 20)  # 0.00017748158796399694
specify_decimal(mean(gamma_025$CCT_includeXsquare), 20)  # 7.927125e-05

specify_decimal(median(gamma_025$SMAT), 20) # 0.00275967613626804109
specify_decimal(median(gamma_025$CCT), 20)  # 4.826216e-06
specify_decimal(median(gamma_025$CCT_includeXsquare), 20)  # 2.721883e-07

specify_decimal(max(gamma_025$SMAT), 20) # 0.78691592406093091583
specify_decimal(max(gamma_025$CCT), 20)  # 0.02158895824613293879
specify_decimal(max(gamma_025$CCT_includeXsquare), 20)  # 0.00913191681141212186

specify_decimal(min(gamma_025$SMAT), 20) # 8.618284e-11
specify_decimal(min(gamma_025$CCT), 20)  # 1.044564e-11
specify_decimal(min(gamma_025$CCT_includeXsquare), 20)  # 8.8818e-16

gamma_035 <- 
  full_df %>%
  filter(Npheno == 12,
         delta == 1,
         delta_select == 0.75,
         gamma == 0.35)
sum(gamma_035$SMAT <= 0.05/50000) # 13
specify_decimal(mean(gamma_035$SMAT), 20)  # 0.02245105569822625907
specify_decimal(mean(gamma_035$CCT), 20)   # 0.00017748158796399339
specify_decimal(mean(gamma_035$CCT_includeXsquare), 20)  # 7.927125e-05, 
# different from gamma = 0.15 or gamma = 0.25 

specify_decimal(median(gamma_035$SMAT), 20) # 0.00275967613626793007
specify_decimal(median(gamma_035$CCT), 20)  # 4.826216e-06
specify_decimal(median(gamma_035$CCT_includeXsquare), 20)  # 2.721883e-07

specify_decimal(max(gamma_035$SMAT), 20) # 0.78691592406092536471
specify_decimal(max(gamma_035$CCT), 20)  # 0.02158895824613282777
specify_decimal(max(gamma_035$CCT_includeXsquare), 20)  # 0.00913191681141212186

specify_decimal(min(gamma_035$SMAT), 20) # 8.618284e-11
specify_decimal(min(gamma_035$CCT), 20)  # 1.044564e-11
specify_decimal(min(gamma_035$CCT_includeXsquare), 20)  # 8.8818e-16
