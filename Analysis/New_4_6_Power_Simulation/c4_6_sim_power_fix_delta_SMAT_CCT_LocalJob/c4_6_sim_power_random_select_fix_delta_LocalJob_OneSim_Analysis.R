# Current array

#' This is revised version of the simulation code
#' This code will be used for simulating the type I error
#' 
# install.packages("MASS")
# install.packages("gdata")
# install.packages("CompQuadForm")
# install.packages("dglm")
# install.packages("lubridate")
# install.packages("gridExtra")
# install.packages("reshape2")

library(gridExtra)
library(reshape2)
library(lubridate)
library(tidyverse)

setwd("~/Dropbox/Emory Courses/DOHG/")
# setwd("/mnt/EpsteinFSS/data/sbian/CauchyTest/")
# setwd("~/Dropbox/Emory Courses/DOHG/")

#################### --------- Laod the data ----------- #######################

###### Predefined global parameter
#' Number of simulations, number of subjects, Number of phenotypes and MAF
#' (independent of covariates)

Nsub_array   <- c(500, 1000, 2000, 5000, 10000)
Npheno_array <- c(4, 6, 8, 10, 12)
maf_array    <- 0.25
gamma_array  <- c(0.15)
delta_array  <- c(0.5)
delta_select_array <- c(1)
Nsim         <- 4

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
  
  file_GAMuT_ProjMatrix_time   <- paste("time", "GAMuT_ProjMatrix_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_GAMuT_LineKernel_time   <- paste("time", "GAMuT_LineKernel_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_SMAT_time               <- paste("time", "SMAT_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_CCT_time                <- paste("time", "CCT_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_includeXsquare_CCT_time <- paste("time", "includeXsquare_CCT_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_simulation_time         <- paste("time", "simulation_time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  
  assign(file_GAMuT_ProjMatrix_time,
         get(load(paste("Result/New_4_6_Power_Simulation_Local_OneSim/", file_GAMuT_ProjMatrix_time, ".RData", sep = ""))))
  
  assign(file_GAMuT_LineKernel_time,
         get(load(paste("Result/New_4_6_Power_Simulation_Local_OneSim/", file_GAMuT_LineKernel_time, ".RData", sep = ""))))

  assign(file_SMAT_time,
         get(load(paste("Result/New_4_6_Power_Simulation_Local_OneSim/", file_SMAT_time, ".RData", sep = ""))))

  assign(file_CCT_time,
         get(load(paste("Result/New_4_6_Power_Simulation_Local_OneSim/", file_CCT_time, ".RData", sep = ""))))

  assign(file_includeXsquare_CCT_time,
         get(load(paste("Result/New_4_6_Power_Simulation_Local_OneSim/", file_includeXsquare_CCT_time, ".RData", sep = ""))))

  assign(file_simulation_time,
         get(load(paste("Result/New_4_6_Power_Simulation_Local_OneSim/", file_simulation_time, ".RData", sep = ""))))
  
  temp_DF <- data.frame(
    "Nsub"          = Nsub,
    "Npheno"        = Npheno,
    "maf"           = maf,
    "gamma"         = gamma,
    "delta"         = delta,
    "delta_select"  = delta_select,
    "GAMuT_ProjMatrix" = get(load(paste("Result/New_4_6_Power_Simulation_Local_OneSim/", file_GAMuT_ProjMatrix_time, ".RData", sep = ""))),
    "GAMuT_LineKernel" = get(load(paste("Result/New_4_6_Power_Simulation_Local_OneSim/", file_GAMuT_LineKernel_time, ".RData", sep = ""))),
    "SMAT"             = get(load(paste("Result/New_4_6_Power_Simulation_Local_OneSim/", file_SMAT_time, ".RData", sep = ""))),
    "CCT"              = get(load(paste("Result/New_4_6_Power_Simulation_Local_OneSim/", file_CCT_time, ".RData", sep = ""))),
    "CCT_includeXsquare" = get(load(paste("Result/New_4_6_Power_Simulation_Local_OneSim/", file_includeXsquare_CCT_time, ".RData", sep = ""))),
    "Overall" = get(load(paste("Result/New_4_6_Power_Simulation_Local_OneSim/", file_simulation_time, ".RData", sep = "")))
  )
  
  full_df <- rbind(full_df, temp_DF)
}

############ ---------- For the Fixed Number of Subjects ----------- ############
# ------ 500 subject
subj500 <- full_df %>%
  filter(Nsub == 500) %>%
  select(Npheno, 
         delta_select,
         GAMuT_ProjMatrix, 
         GAMuT_LineKernel,
         SMAT,
         CCT,
         CCT_includeXsquare) %>%
  gather(., "Method", "Value", GAMuT_ProjMatrix:CCT_includeXsquare) %>% 
  mutate(NphenoF = factor(Npheno, 
                          levels = c(4, 6, 8, 10, 12), 
                          labels = c("4", "6", "8", "10", "12"))) %>%
  mutate(MethodF = factor(Method, 
                         levels = c("CCT", "CCT_includeXsquare", "SMAT", "GAMuT_LineKernel", "GAMuT_ProjMatrix"), 
                          labels = c("CCT", "CCT (with X^2)", "SMAT", "GAMuT LineKernel", "GAMuT ProjMatrix"))) %>%
  select(-Npheno, -Method) %>%
  group_by(MethodF, NphenoF) %>%
  summarise(mean = mean(Value))

g <- ggplot(subj500, aes(x = NphenoF, y = mean, color = MethodF, group = MethodF)) + 
  geom_point() +
  geom_path(aes(colour = MethodF, group = MethodF)) +
  ggtitle("Time per Power Simulation for 500 Subjects") +
  xlab("Number of phenotypes") + 
  ylab("Time in Second") + 
  theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  ) +
  ylim(0, 40) +
  theme_bw() +
  theme(text = element_text(size = 20)) 
gt <- tableGrob(subj500 %>% 
                  mutate(mean = round(mean, 3)) %>% 
                  spread(., NphenoF, mean), 
                rows = NULL, theme = ttheme_minimal())
png("Analysis/New_4_6_Power_Simulation/New_4_6_Power_Simulation_Local_OneSim/Given500Subj_v2.png",
    width = 1000, height = 600)
grid.arrange(g, gt, ncol = 1, heights=c(2,1))
dev.off()

# ------ 1000 subject
subj1000 <- full_df %>%
  filter(Nsub == 1000) %>%
  select(Npheno, 
         delta_select,
         GAMuT_ProjMatrix, 
         GAMuT_LineKernel,
         SMAT,
         CCT,
         CCT_includeXsquare) %>%
  gather(., "Method", "Value", GAMuT_ProjMatrix:CCT_includeXsquare) %>% 
  mutate(NphenoF = factor(Npheno, 
                          levels = c(4, 6, 8, 10, 12), 
                          labels = c("4", "6", "8", "10", "12"))) %>%
  mutate(MethodF = factor(Method, 
                          levels = c("CCT", "CCT_includeXsquare", "SMAT", "GAMuT_LineKernel", "GAMuT_ProjMatrix"), 
                          labels = c("CCT", "CCT (with X^2)", "SMAT", "GAMuT LineKernel", "GAMuT ProjMatrix"))) %>%
  select(-Npheno, -Method) %>%
  group_by(MethodF, NphenoF) %>%
  summarise(mean = mean(Value))

g <- ggplot(subj1000, aes(x = NphenoF, y = mean, color = MethodF, group = MethodF)) + 
  geom_point() +
  geom_path(aes(colour = MethodF, group = MethodF)) +
  ggtitle("Time per Power Simulation for 1,000 Subjects") +
  xlab("Number of phenotypes") + 
  ylab("Time in Second") + 
  theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  ) +
  ylim(0, 40) +
  theme_bw() +
  theme(text = element_text(size = 20)) 
gt <- tableGrob(subj1000 %>% 
                  mutate(mean = round(mean, 3)) %>% 
                  spread(., NphenoF, mean), 
                rows = NULL, theme = ttheme_minimal())
png("Analysis/New_4_6_Power_Simulation/New_4_6_Power_Simulation_Local_OneSim/Given1000Subj_v2.png",
    width = 1000, height = 600)
grid.arrange(g, gt, ncol = 1, heights=c(2,1))
dev.off()

# ------ 2000 subject
subj2000 <- full_df %>%
  filter(Nsub == 2000) %>%
  select(Npheno, 
         delta_select,
         GAMuT_ProjMatrix, 
         GAMuT_LineKernel,
         SMAT,
         CCT,
         CCT_includeXsquare) %>%
  gather(., "Method", "Value", GAMuT_ProjMatrix:CCT_includeXsquare) %>% 
  mutate(NphenoF = factor(Npheno, 
                          levels = c(4, 6, 8, 10, 12), 
                          labels = c("4", "6", "8", "10", "12"))) %>%
  mutate(MethodF = factor(Method, 
                          levels = c("CCT", "CCT_includeXsquare", "SMAT", "GAMuT_LineKernel", "GAMuT_ProjMatrix"), 
                          labels = c("CCT", "CCT (with X^2)", "SMAT", "GAMuT LineKernel", "GAMuT ProjMatrix"))) %>%
  select(-Npheno, -Method) %>%
  group_by(MethodF, NphenoF) %>%
  summarise(mean = mean(Value))

g <- ggplot(subj2000, aes(x = NphenoF, y = mean, color = MethodF, group = MethodF)) + 
  geom_point() +
  geom_path(aes(colour = MethodF, group = MethodF)) +
  ggtitle("Time per Power Simulation for 2,000 Subjects") +
  xlab("Number of phenotypes") + 
  ylab("Time in Second") + 
  theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  ) +
  ylim(0, 40) +
  theme_bw()+
  theme(text = element_text(size = 20)) 
gt <- tableGrob(subj2000 %>% 
                  mutate(mean = round(mean, 3)) %>% 
                  spread(., NphenoF, mean), 
                rows = NULL, theme = ttheme_minimal())
png("Analysis/New_4_6_Power_Simulation/New_4_6_Power_Simulation_Local_OneSim/Given2000Subj_v2.png",
    width = 1000, height = 600)
grid.arrange(g, gt, ncol = 1, heights=c(2,1))
dev.off()

# ------ 5000 subject
subj5000 <- full_df %>%
  filter(Nsub == 5000) %>%
  select(Npheno, 
         delta_select,
         GAMuT_ProjMatrix, 
         GAMuT_LineKernel,
         SMAT,
         CCT,
         CCT_includeXsquare) %>%
  gather(., "Method", "Value", GAMuT_ProjMatrix:CCT_includeXsquare) %>% 
  mutate(NphenoF = factor(Npheno, 
                          levels = c(4, 6, 8, 10, 12), 
                          labels = c("4", "6", "8", "10", "12"))) %>%
  mutate(MethodF = factor(Method, 
                          levels = c("CCT", "CCT_includeXsquare", "SMAT", "GAMuT_LineKernel", "GAMuT_ProjMatrix"), 
                          labels = c("CCT", "CCT (with X^2)", "SMAT", "GAMuT LineKernel", "GAMuT ProjMatrix"))) %>%
  select(-Npheno, -Method) %>%
  group_by(MethodF, NphenoF) %>%
  summarise(mean = mean(Value))

g <- ggplot(subj5000, aes(x = NphenoF, y = mean, color = MethodF, group = MethodF)) + 
  geom_point() +
  geom_path(aes(colour = MethodF, group = MethodF)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, seq(1, 16, by = 2), seq(18, 150, by = 50),  seq(200, 550, 50)),
                     minor_breaks = seq(18, 150, by = 50)) +
  ggtitle("Time per Power Simulation for 5,000 Subjects") +
  xlab("Number of phenotypes") + 
  ylab("Time in Second") + 
  theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  ) +
  theme_bw()+
  theme(text = element_text(size = 20)) 
gt <- tableGrob(subj5000 %>% 
                  mutate(mean = round(mean, 3)) %>% 
                  spread(., NphenoF, mean), 
                rows = NULL, theme = ttheme_minimal())
png("Analysis/New_4_6_Power_Simulation/New_4_6_Power_Simulation_Local_OneSim/Given5000Subj_v2.png",
    width = 1000, height = 600)
grid.arrange(g, gt, ncol = 1, heights=c(2,1))
dev.off()

# ------ 10000 subject
subj10000 <- full_df %>%
  filter(Nsub == 10000) %>%
  select(Npheno, 
         delta_select,
         GAMuT_ProjMatrix, 
         GAMuT_LineKernel,
         SMAT,
         CCT,
         CCT_includeXsquare) %>%
  gather(., "Method", "Value", GAMuT_ProjMatrix:CCT_includeXsquare) %>% 
  mutate(NphenoF = factor(Npheno, 
                          levels = c(4, 6, 8, 10, 12), 
                          labels = c("4", "6", "8", "10", "12"))) %>%
  mutate(MethodF = factor(Method, 
                          levels = c("CCT", "CCT_includeXsquare", "SMAT", "GAMuT_LineKernel", "GAMuT_ProjMatrix"), 
                          labels = c("CCT", "CCT (with X^2)", "SMAT", "GAMuT LineKernel", "GAMuT ProjMatrix"))) %>%
  select(-Npheno, -Method) %>%
  group_by(MethodF, NphenoF) %>%
  summarise(mean = mean(Value))

g <- ggplot(subj10000, aes(x = NphenoF, y = mean, color = MethodF, group = MethodF)) + 
  geom_point() +
  geom_path(aes(colour = MethodF, group = MethodF)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, seq(1, 32, by = 2), seq(34, 2100, by = 500),  seq(2150, 2300, 50)),
                     minor_breaks = seq(18, 150, by = 50)) +
  ggtitle("Time per Power Simulation for 10,000 Subjects") +
  xlab("Number of phenotypes") + 
  ylab("Time in Second") + 
  theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  ) +
  theme_bw()+
  theme(text = element_text(size = 20)) 
gt <- tableGrob(subj10000 %>% 
                  mutate(mean = round(mean, 3)) %>% 
                  spread(., NphenoF, mean), 
                rows = NULL, theme = ttheme_minimal())
png("Analysis/New_4_6_Power_Simulation/New_4_6_Power_Simulation_Local_OneSim/Given10000Subj_v2.png",
    width = 1000, height = 600)
grid.arrange(g, gt, ncol = 1, heights=c(2,1))
dev.off()




########### ---------- For the Fixed Number of Phenotypes ----------- ##########
# ------ 4 phenotypes
pheno4 <- full_df %>%
  filter(Npheno == 4) %>%
  select(Nsub, 
         delta_select,
         GAMuT_ProjMatrix, 
         GAMuT_LineKernel,
         SMAT,
         CCT,
         CCT_includeXsquare) %>%
  gather(., "Method", "Value", GAMuT_ProjMatrix:CCT_includeXsquare) %>% 
  mutate(MethodF = factor(Method, 
                          levels = c("CCT", "CCT_includeXsquare", "SMAT", "GAMuT_LineKernel", "GAMuT_ProjMatrix"), 
                          labels = c("CCT", "CCT (with X^2)", "SMAT", "GAMuT LineKernel", "GAMuT ProjMatrix"))) %>%
  select(-Method) %>%
  filter(!MethodF %in% c("GAMuT LineKernel", "GAMuT ProjMatrix")) %>%
  group_by(MethodF, Nsub) %>%
  summarise(mean = mean(Value), .groups = 'drop')
g <- ggplot(pheno4, aes(x = Nsub, y = mean, color = MethodF, group = MethodF)) + 
  geom_point() +
  geom_path(aes(colour = MethodF, group = MethodF)) +
  ggtitle("Time per Power Simulation for 4 Phenotypes") +
  xlab("Number of subjects") + 
  ylab("Time in Second") + 
  theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  ) +
  ylim(0, 40) +
  scale_x_continuous(breaks = c(500, 1000, 2000, 5000, 10000),
                     labels = c("500", "1,000", "2,000", "5,000", "10,000")) +
  theme_bw() +
  theme(text = element_text(size = 20)) 
gt <- tableGrob(pheno4 %>% 
                  mutate(mean = round(mean, 3)) %>% 
                  spread(., Nsub, mean), 
                rows = NULL, theme = ttheme_minimal())
png("Analysis/New_4_6_Power_Simulation/New_4_6_Power_Simulation_Local_OneSim/Given4Pheno_v2.png",
    width = 1000, height = 600)
grid.arrange(g, gt, ncol = 1, heights=c(2,1))
dev.off()

# ------ 6 phenotypes
pheno6 <- full_df %>%
  filter(Npheno == 6) %>%
  select(Nsub, 
         delta_select,
         GAMuT_ProjMatrix, 
         GAMuT_LineKernel,
         SMAT,
         CCT,
         CCT_includeXsquare) %>%
  gather(., "Method", "Value", GAMuT_ProjMatrix:CCT_includeXsquare) %>% 
  mutate(MethodF = factor(Method, 
                          levels = c("CCT", "CCT_includeXsquare", "SMAT", "GAMuT_LineKernel", "GAMuT_ProjMatrix"), 
                          labels = c("CCT", "CCT (with X^2)", "SMAT", "GAMuT LineKernel", "GAMuT ProjMatrix"))) %>%
  select(-Method) %>%
  filter(!MethodF %in% c("GAMuT LineKernel", "GAMuT ProjMatrix")) %>%
  group_by(MethodF, Nsub) %>%
  summarise(mean = mean(Value), .groups = 'drop')

g <- ggplot(pheno6, aes(x = Nsub, y = mean, color = MethodF, group = MethodF)) + 
  geom_point() +
  geom_path(aes(colour = MethodF, group = MethodF)) +
  ggtitle("Time per Power Simulation for 6 Phenotypes") +
  xlab("Number of subjects") + 
  ylab("Time in Second") + 
  theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  ) +
  ylim(0, 40) +
  scale_x_continuous(breaks = c(500, 1000, 2000, 5000, 10000),
                     labels = c("500", "1,000", "2,000", "5,000", "10,000")) +
  theme_bw() +
  theme(text = element_text(size = 20)) 
gt <- tableGrob(pheno6 %>% 
                  mutate(mean = round(mean, 3)) %>% 
                  spread(., Nsub, mean), 
                rows = NULL, theme = ttheme_minimal())
png("Analysis/New_4_6_Power_Simulation/New_4_6_Power_Simulation_Local_OneSim/Given6Pheno_v2.png",
    width = 1000, height = 600)
grid.arrange(g, gt, ncol = 1, heights=c(2,1))
dev.off()

# ------ 8 phenotypes
pheno8 <- full_df %>%
  filter(Npheno == 8) %>%
  select(Nsub, 
         delta_select,
         GAMuT_ProjMatrix, 
         GAMuT_LineKernel,
         SMAT,
         CCT,
         CCT_includeXsquare) %>%
  gather(., "Method", "Value", GAMuT_ProjMatrix:CCT_includeXsquare) %>% 
  mutate(MethodF = factor(Method, 
                          levels = c("CCT", "CCT_includeXsquare", "SMAT", "GAMuT_LineKernel", "GAMuT_ProjMatrix"), 
                          labels = c("CCT", "CCT (with X^2)", "SMAT", "GAMuT LineKernel", "GAMuT ProjMatrix"))) %>%
  select(-Method) %>%
  filter(!MethodF %in% c("GAMuT LineKernel", "GAMuT ProjMatrix")) %>%
  group_by(MethodF, Nsub) %>%
  summarise(mean = mean(Value), .groups = 'drop')
g <- ggplot(pheno8, aes(x = Nsub, y = mean, color = MethodF, group = MethodF)) + 
  geom_point() +
  geom_path(aes(colour = MethodF, group = MethodF)) +
  ggtitle("Time per Power Simulation for 8 Phenotypes") +
  xlab("Number of subjects") + 
  ylab("Time in Second") + 
  theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  ) +
  ylim(0, 40) +
  scale_x_continuous(breaks = c(500, 1000, 2000, 5000, 10000),
                     labels = c("500", "1,000", "2,000", "5,000", "10,000")) +
  theme_bw() +
  theme(text = element_text(size = 20)) 
gt <- tableGrob(pheno8 %>% 
                  mutate(mean = round(mean, 3)) %>% 
                  spread(., Nsub, mean), 
                rows = NULL, theme = ttheme_minimal())
png("Analysis/New_4_6_Power_Simulation/New_4_6_Power_Simulation_Local_OneSim/Given8Pheno_v2.png",
    width = 1000, height = 600)
grid.arrange(g, gt, ncol = 1, heights=c(2,1))
dev.off()

# ------ 10 phenotypes
pheno10 <- full_df %>%
  filter(Npheno == 10) %>%
  select(Nsub, 
         delta_select,
         GAMuT_ProjMatrix, 
         GAMuT_LineKernel,
         SMAT,
         CCT,
         CCT_includeXsquare) %>%
  gather(., "Method", "Value", GAMuT_ProjMatrix:CCT_includeXsquare) %>% 
  mutate(MethodF = factor(Method, 
                          levels = c("CCT", "CCT_includeXsquare", "SMAT", "GAMuT_LineKernel", "GAMuT_ProjMatrix"), 
                          labels = c("CCT", "CCT (with X^2)", "SMAT", "GAMuT LineKernel", "GAMuT ProjMatrix"))) %>%
  select(-Method) %>%
  filter(!MethodF %in% c("GAMuT LineKernel", "GAMuT ProjMatrix")) %>%
  ungroup() %>%
  select(Nsub, MethodF, Value) %>%
  group_by(MethodF, Nsub) %>%
  summarise(mean = mean(Value), .groups = 'drop')

g <- ggplot(pheno10, aes(x = Nsub, y = mean, color = MethodF, group = MethodF)) + 
  geom_point() +
  geom_path(aes(colour = MethodF, group = MethodF)) +
  ggtitle("Time per Power Simulation for 10 Phenotypes") +
  xlab("Number of subjects") + 
  ylab("Time in Second") + 
  theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  ) +
  ylim(0, 40) +
  scale_x_continuous(breaks = c(500, 1000, 2000, 5000, 10000),
                     labels = c("500", "1,000", "2,000", "5,000", "10,000")) +
  theme_bw() +
  theme(text = element_text(size = 20)) 
gt <- tableGrob(pheno10 %>% 
                  mutate(mean = round(mean, 3)) %>% 
                  spread(., Nsub, mean), 
                rows = NULL, theme = ttheme_minimal())
png("Analysis/New_4_6_Power_Simulation/New_4_6_Power_Simulation_Local_OneSim/Given10Pheno_v2.png",
    width = 1000, height = 600)
grid.arrange(g, gt, ncol = 1, heights=c(2,1))
dev.off()

# ------ 12 phenotypes
pheno12 <- full_df %>%
  filter(Npheno == 12) %>%
  select(Nsub, 
         delta_select,
         GAMuT_ProjMatrix, 
         GAMuT_LineKernel,
         SMAT,
         CCT,
         CCT_includeXsquare) %>%
  gather(., "Method", "Value", GAMuT_ProjMatrix:CCT_includeXsquare) %>% 
  mutate(MethodF = factor(Method, 
                          levels = c("CCT", "CCT_includeXsquare", "SMAT", "GAMuT_LineKernel", "GAMuT_ProjMatrix"), 
                          labels = c("CCT", "CCT (with X^2)", "SMAT", "GAMuT LineKernel", "GAMuT ProjMatrix"))) %>%
  select(-Method) %>%
  filter(!MethodF %in% c("GAMuT LineKernel", "GAMuT ProjMatrix")) %>%
  select(Nsub, MethodF, Value) %>%
  group_by(MethodF, Nsub) %>%
  summarise(mean = mean(Value), .groups = 'drop')

g <- ggplot(pheno12, aes(x = Nsub, y = mean, color = MethodF, group = MethodF)) + 
  geom_point() +
  geom_path(aes(colour = MethodF, group = MethodF)) +
  ggtitle("Time per Power Simulation for 12 Phenotypes") +
  xlab("Number of subjects") + 
  ylab("Time in Second") + 
  theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  ) +
  ylim(0, 40) +
  scale_x_continuous(breaks = c(500, 1000, 2000, 5000, 10000),
                     labels = c("500", "1,000", "2,000", "5,000", "10,000")) +
  theme_bw() +
  theme(text = element_text(size = 20)) 
gt <- tableGrob(pheno12 %>% 
                  mutate(mean = round(mean, 3)) %>% 
                  spread(., Nsub, mean), 
                rows = NULL, theme = ttheme_minimal())
png("Analysis/New_4_6_Power_Simulation/New_4_6_Power_Simulation_Local_OneSim/Given12Pheno_v2.png",
    width = 1000, height = 600)
grid.arrange(g, gt, ncol = 1, heights=c(2,1))
dev.off()

# Bullet Point 1
2.590*100000/3600/24
2.298*100000/3600/24
30.094*100000/3600/24

# Simple linear regression
datSMAT <- pheno12[pheno12$MethodF=="SMAT",]
SMAT_LM <- lm(mean~Nsub, data = datSMAT)
0.40564 + 0.00302 * 50000

datCCT <- pheno12[pheno12$MethodF=="CCT",]
CCT_LM <- lm(mean~Nsub, data = datCCT)
0.0727902 + 0.0002512 * 50000

datCCT_X2 <- pheno12[pheno12$MethodF=="CCT (with X^2)",]
CCT_X2_LM <- lm(mean~Nsub, data = datCCT_X2)
0.1415871 + 0.0002162 * 50000

10.95159*100000/3600/24
