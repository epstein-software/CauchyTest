library(tidyverse)

setwd("~/Dropbox/Emory Courses/DOHG/")


###### Predefined global parameter
#' Number of simulations, number of subjects, Number of phenotypes and MAF
#' (independent of covariates)

Nsub_array   <- c(5000, 10000)
Npheno_array <- c(4, 6, 8, 10, 12)
maf_array    <- 0.25
gamma_array  <- c(NA, 0.15, 0.55, 0.95)
delta_array  <- c(0.35, 0.45, 0.6, 1.5)
delta_select_array <- c(0.25, 0.5, 0.75, 1)
Nsim         <- 10000

param_comb <- expand.grid(Nsub   = Nsub_array,
                          Npheno = Npheno_array,
                          maf    = maf_array,
                          gamma  = gamma_array,
                          delta  = delta_array,
                          delta_select = delta_select_array)

############ ---------- Construct the Overall Dataframe ----------- ############

full_pvalue_df <- data.frame()
for (array.id in 1:nrow(param_comb)) {
  
  Nsub   <- param_comb[array.id, "Nsub"]
  Npheno <- param_comb[array.id, "Npheno"]
  maf    <- param_comb[array.id, "maf"]
  gamma  <- param_comb[array.id, "gamma"]
  delta  <- param_comb[array.id, "delta"]
  delta_select  <- param_comb[array.id, "delta_select"]
  
  file_name_SMAT_pvalue               <- paste("time_power", "SMAT", "time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_name_includeXsquare_SMAT_pvalue <- paste("time_power", "includeXsquare_SMAT", "time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_name_CCT_pvalue                <- paste("time_power", "CCT", "time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_name_includeXsquare_CCT_pvalue <- paste("time_power", "includeXsquare_CCT", "time", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  
  assign(file_name_SMAT_pvalue,
         get(load(paste("Result/New_4_6_Power_Simulation_ExcludeGAMUT/", file_name_SMAT_pvalue, ".RData", sep = ""))))
  
  assign(file_name_includeXsquare_SMAT_pvalue,
         get(load(paste("Result/New_4_6_Power_Simulation_ExcludeGAMUT/", file_name_includeXsquare_SMAT_pvalue, ".RData", sep = ""))))
  
  assign(file_name_CCT_pvalue,
         get(load(paste("Result/New_4_6_Power_Simulation_ExcludeGAMUT/", file_name_CCT_pvalue, ".RData", sep = ""))))
  
  assign(file_name_includeXsquare_CCT_pvalue,
         get(load(paste("Result/New_4_6_Power_Simulation_ExcludeGAMUT/", file_name_includeXsquare_CCT_pvalue, ".RData", sep = ""))))
  
  temp_DF <- data.frame(
    "Nsub"          = Nsub,
    "Npheno"        = Npheno,
    "maf"           = maf,
    "gamma"         = gamma,
    "delta"         = delta,
    "delta_select"  = delta_select,
    "SMAT"                = get(load(paste("Result/New_4_6_Power_Simulation_ExcludeGAMUT/", file_name_SMAT_pvalue, ".RData", sep = ""))),
    "SMAT_includeXsquare" = get(load(paste("Result/New_4_6_Power_Simulation_ExcludeGAMUT/", file_name_includeXsquare_SMAT_pvalue, ".RData", sep = ""))),
    "CCT"                 = get(load(paste("Result/New_4_6_Power_Simulation_ExcludeGAMUT/", file_name_CCT_pvalue, ".RData", sep = ""))),
    "CCT_includeXsquare"  = get(load(paste("Result/New_4_6_Power_Simulation_ExcludeGAMUT/", file_name_includeXsquare_CCT_pvalue, ".RData", sep = "")))
  )
  
  full_pvalue_df <- rbind(full_pvalue_df, temp_DF)
}

############### Fix number of subject, gamma and delta #########################
temp_df <- full_pvalue_df %>%
  select(-maf) 

temp_df_F <- temp_df %>%
  select(-gamma, -delta, -delta_select) %>% # 6400000       6
  gather(., 
         key = "Method", 
         value = "Power",
         SMAT,
         SMAT_includeXsquare,
         CCT,
         CCT_includeXsquare) %>%
  group_by(., 
           Nsub,
           Npheno,
           Method) %>%
  summarise(Ave_Time = mean(Power)) 
temp_df_F$Npheno_F = factor(temp_df_F$Npheno)
table(temp_df_F$Npheno, temp_df_F$Npheno_F)

ggplot(temp_df_F, 
       aes(x = Npheno_F, y = Ave_Time,  group=Method)) +
  geom_line(aes(linetype=Method, color = Method))+
  theme_bw() +
  facet_wrap(~ Nsub, ncol = 2) + 
  theme(strip.text = element_text(size=12, face="bold", color="darkblue")) +
  ggtitle("Average time for completing one single simulation: an average from 6,400,000 simulations") +
  xlab("Number of phenotypes") + 
  ylab("Simulation time in second") +
  theme(
    plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=16, face="bold"),
    axis.title.y = element_text(color="#993333", size=16, face="bold")
  ) 
ggsave(paste("Analysis/Power_Time_Simulation_Dec2021/Simulation_Time_01172022.pdf", sep = ""),
       width = 40, height = 20, units = "cm")

temp_df_F$`Time for 10,000 Simulations in Hour` = temp_df_F$Ave_Time*10000/3600
names(temp_df_F)[names(temp_df_F)=="Ave_Time"] = "Average time per simulation in second"
write_csv(temp_df_F %>% select(-Npheno_F), 
          "Analysis/Power_Time_Simulation_Dec2021/Simulation_Time_01182022.csv")

########### Only without sparsity
temp_df <- full_pvalue_df %>%
  filter(delta_select == 1) %>%
  select(-maf) 

temp_df_F <- temp_df %>%
  select(-gamma, -delta, -delta_select) %>% # 1600000       6
  gather(., 
         key = "Method", 
         value = "Power",
         SMAT,
         SMAT_includeXsquare,
         CCT,
         CCT_includeXsquare) %>%
  group_by(., 
           Nsub,
           Npheno,
           Method) %>%
  summarise(Ave_Time = mean(Power)) 
temp_df_F$Npheno_F = factor(temp_df_F$Npheno)
table(temp_df_F$Npheno, temp_df_F$Npheno_F)

ggplot(temp_df_F, 
       aes(x = Npheno_F, y = Ave_Time,  group=Method)) +
  geom_line(aes(linetype=Method, color = Method))+
  theme_bw() +
  facet_wrap(~ Nsub, ncol = 2) + 
  theme(strip.text = element_text(size=12, face="bold", color="darkblue")) +
  ggtitle("Average time for completing one single simulation: an average from 1,600,000 simulations, only include the simulation without sparsity") +
  xlab("Number of phenotypes") + 
  ylab("Simulation time in second") +
  theme(
    plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=16, face="bold"),
    axis.title.y = element_text(color="#993333", size=16, face="bold")
  ) 
ggsave(paste("Analysis/Power_Time_Simulation_Dec2021/Simulation_Time_wo_sparsity_01172022.pdf", sep = ""),
       width = 40, height = 20, units = "cm")

temp_df_F$`Time for 10,000 Simulations in Hour` = temp_df_F$Ave_Time*10000/3600
names(temp_df_F)[names(temp_df_F)=="Ave_Time"] = "Average time per simulation in second"
write_csv(temp_df_F %>% select(-Npheno_F), 
          "Analysis/Power_Time_Simulation_Dec2021/Simulation_Time_wo_sparsity_01172022.csv")

