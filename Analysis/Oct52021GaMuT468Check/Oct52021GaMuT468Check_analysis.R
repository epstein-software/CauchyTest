library(tidyverse)
setwd("~/Dropbox/Emory Courses/DOHG")
source("Analysis/result_analysis_source_code.R")

############################# 2,000 Simulation #############################

## ---- 500 subjects + 2000 simulations ---- ###
all_files <- list.files("Result/Oct52021GaMuT468Check/")
subject_500_dglm_2000 <- grep("*._GAMuT_time_adj_500_subject_dglm_2000.*", all_files, value = T)

for (i in 1:length(subject_500_dglm_2000)) {
  file_name <- subject_500_dglm_2000[i]
  load(paste("Result/Oct52021GaMuT468Check/", file_name, sep = ""))
  pheno = substr(file_name, nchar(file_name) - 7, nchar(file_name) - 6)
  assign(paste("subject_500_dglm_2000", pheno, "_Pheno", sep = ""), GAMuT_time)
}

CCT_pvalue_df_1 <- data.frame(Number_Sub = 500,
                            Number_Sim = 2000,
                            subject_500_dglm_2000_4_Pheno,
                            subject_500_dglm_2000_6_Pheno,
                            subject_500_dglm_2000_8_Pheno)
apply(CCT_pvalue_df_1[, 3:5], 2, mean)
# subject_500_dglm_2000_4_Pheno subject_500_dglm_2000_6_Pheno subject_500_dglm_2000_8_Pheno 
#                     0.4808831                     0.5027703                     0.5221191 
names(CCT_pvalue_df_1)[3:5] <- c("4_pheno", "6_pheno", "8_pheno")

## ---- 1000 subjects + 2000 simulations ---- ###
all_files <- list.files("Result/Oct52021GaMuT468Check/")
subject_1000_dglm_2000 <- grep("*._GAMuT_time_adj_1000_subject_dglm_2000.*", all_files, value = T)

for (i in 1:length(subject_1000_dglm_2000)) {
  file_name <- subject_1000_dglm_2000[i]
  load(paste("Result/Oct52021GaMuT468Check/", file_name, sep = ""))
  pheno = substr(file_name, nchar(file_name) - 7, nchar(file_name) - 6)
  assign(paste("subject_1000_dglm_2000", pheno, "_Pheno", sep = ""), GAMuT_time)
}

CCT_pvalue_df_2 <- data.frame(Number_Sub = 1000,
                            Number_Sim = 2000,
                            subject_1000_dglm_2000_4_Pheno,
                            subject_1000_dglm_2000_6_Pheno,
                            subject_1000_dglm_2000_8_Pheno)
apply(CCT_pvalue_df_2[, 3:5], 2, mean)
# subject_1000_dglm_2000_4_Pheno subject_1000_dglm_2000_6_Pheno subject_1000_dglm_2000_8_Pheno 
#                       3.390311                       3.421680                       3.660177 
names(CCT_pvalue_df_2)[3:5] <- c("4_pheno", "6_pheno", "8_pheno")

df_2000_simu <- rbind(CCT_pvalue_df_1, 
                      CCT_pvalue_df_2) 

temp = 
  df_2000_simu %>%
  dplyr::select(., -Number_Sim) %>%
  gather(., "Number Phenotypes", "Computation Time", -Number_Sub)%>%
  mutate(Number_SubF = factor(Number_Sub))
table(temp$Number_Sub, temp$Number_SubF)

# ---- Median
ggsave("Analysis/Oct52021GaMuT468Check/median_computationalTime_over_numberPheno_2000simu.pdf", width = 12, height = 6)
df_2000_simu %>%
  dplyr::select(., -Number_Sim) %>%
  gather(., "Number Phenotypes", "Computation Time", - Number_Sub)%>%
  mutate("Number Subjects" = factor(Number_Sub)) %>%
  dplyr::select(., -Number_Sub) %>%
  group_by(`Number Subjects`, `Number Phenotypes`) %>%
  summarise(Median = median(`Computation Time`)) %>%
  ggplot(., aes(x=`Number Phenotypes`, y=Median, fill= `Number Subjects`)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label = round(Median,2)), size = 3,  position = position_dodge(0.9)) +
  scale_fill_brewer(direction = 1) +
  theme_minimal() +
  labs(x = "Number of Phenotypes",
       y = "Median of Computation Time (Second)",
       title="Median of Computational Time for Different Number of Subjects and Phenotypes",
       subtitle = "2,000 simulations are applied") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10)) +
  theme(plot.title = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) 
dev.off()

# ---- Mean
ggsave("Analysis//Oct52021GaMuT468Check/mean_computationalTime_over_numberPheno_2000simu.pdf", width = 12, height = 6)
df_2000_simu %>%
  dplyr::select(., -Number_Sim) %>%
  gather(., "Number Phenotypes", "Computation Time", - Number_Sub)%>%
  mutate("Number Subjects" = factor(Number_Sub)) %>%
  dplyr::select(., -Number_Sub) %>%
  group_by(`Number Subjects`, `Number Phenotypes`) %>%
  summarise(Mean = mean(`Computation Time`)) %>%
  ggplot(., aes(x=`Number Phenotypes`, y=Mean, fill= `Number Subjects`)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label = round(Mean,2)), size = 3,  position = position_dodge(0.9)) +
  scale_fill_brewer(direction = 1) +
  theme_minimal() +
  labs(x = "Number of Phenotypes",
       y = "Mean of Computation Time (Second)",
       title="Mean of Computational Time for Different Number of Subjects and Phenotypes",
       subtitle = "2,000 simulations are applied") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10)) +
  theme(plot.title = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) 
dev.off()


############################# 1,000 Simulation #############################
rm(list = ls())
## ---- 500 subjects + 1000 simulations ---- ###
all_files <- list.files("Result/Oct52021GaMuT468Check/")
subject_500_dglm_2000 <- grep("*._GAMuT_time_adj_500_subject_dglm_1000.*", all_files, value = T)

for (i in 1:length(subject_500_dglm_2000)) {
  file_name <- subject_500_dglm_2000[i]
  load(paste("Result/Oct52021GaMuT468Check/", file_name, sep = ""))
  pheno = substr(file_name, nchar(file_name) - 7, nchar(file_name) - 6)
  assign(paste("subject_500_dglm_2000", pheno, "_Pheno", sep = ""), GAMuT_time)
}

CCT_pvalue_df_1 <- data.frame(Number_Sub = 500,
                              Number_Sim = 1000,
                              subject_500_dglm_2000_4_Pheno,
                              subject_500_dglm_2000_6_Pheno,
                              subject_500_dglm_2000_8_Pheno)
apply(CCT_pvalue_df_1[, 3:5], 2, mean)
# subject_500_dglm_2000_4_Pheno subject_500_dglm_2000_6_Pheno subject_500_dglm_2000_8_Pheno 
#                    0.4810730                     0.5103602                     0.5323802 
names(CCT_pvalue_df_1)[3:5] <- c("4_pheno", "6_pheno", "8_pheno")

## ---- 1000 subjects + 2000 simulations ---- ###
all_files <- list.files("Result/Oct52021GaMuT468Check/")
subject_1000_dglm_2000 <- grep("*._GAMuT_time_adj_1000_subject_dglm_1000.*", all_files, value = T)

for (i in 1:length(subject_1000_dglm_2000)) {
  file_name <- subject_1000_dglm_2000[i]
  load(paste("Result/Oct52021GaMuT468Check/", file_name, sep = ""))
  pheno = substr(file_name, nchar(file_name) - 7, nchar(file_name) - 6)
  assign(paste("subject_1000_dglm_2000", pheno, "_Pheno", sep = ""), GAMuT_time)
}

CCT_pvalue_df_2 <- data.frame(Number_Sub = 1000,
                              Number_Sim = 1000,
                              subject_1000_dglm_2000_4_Pheno,
                              subject_1000_dglm_2000_6_Pheno,
                              subject_1000_dglm_2000_8_Pheno)
apply(CCT_pvalue_df_2[, 3:5], 2, mean)
# subject_1000_dglm_2000_4_Pheno subject_1000_dglm_2000_6_Pheno subject_1000_dglm_2000_8_Pheno 
#                     3.455465                       3.483398                       3.520599 
names(CCT_pvalue_df_2)[3:5] <- c("4_pheno", "6_pheno", "8_pheno")

df_2000_simu <- rbind(CCT_pvalue_df_1, 
                      CCT_pvalue_df_2) 

temp = 
  df_2000_simu %>%
  dplyr::select(., -Number_Sim) %>%
  gather(., "Number Phenotypes", "Computation Time", -Number_Sub)%>%
  mutate(Number_SubF = factor(Number_Sub))
table(temp$Number_Sub, temp$Number_SubF)

# ---- Median
ggsave("Analysis//Oct52021GaMuT468Check/median_computationalTime_over_numberPheno_1000simu.pdf", width = 12, height = 6)
df_2000_simu %>%
  dplyr::select(., -Number_Sim) %>%
  gather(., "Number Phenotypes", "Computation Time", - Number_Sub)%>%
  mutate("Number Subjects" = factor(Number_Sub)) %>%
  dplyr::select(., -Number_Sub) %>%
  group_by(`Number Subjects`, `Number Phenotypes`) %>%
  summarise(Median = median(`Computation Time`)) %>%
  ggplot(., aes(x=`Number Phenotypes`, y=Median, fill= `Number Subjects`)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label = round(Median,2)), size = 3,  position = position_dodge(0.9)) +
  scale_fill_brewer(direction = 1) +
  theme_minimal() +
  labs(x = "Number of Phenotypes",
       y = "Median of Computation Time (Second)",
       title="Median of Computational Time for Different Number of Subjects and Phenotypes",
       subtitle = "1,000 simulations are applied") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10)) +
  theme(plot.title = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) 
dev.off()

# ---- Mean
ggsave("Analysis//Oct52021GaMuT468Check/mean_computationalTime_over_numberPheno_1000simu.pdf", width = 12, height = 6)
df_2000_simu %>%
  dplyr::select(., -Number_Sim) %>%
  gather(., "Number Phenotypes", "Computation Time", - Number_Sub)%>%
  mutate("Number Subjects" = factor(Number_Sub)) %>%
  dplyr::select(., -Number_Sub) %>%
  group_by(`Number Subjects`, `Number Phenotypes`) %>%
  summarise(Mean = mean(`Computation Time`)) %>%
  ggplot(., aes(x=`Number Phenotypes`, y=Mean, fill= `Number Subjects`)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label = round(Mean,2)), size = 3,  position = position_dodge(0.9)) +
  scale_fill_brewer(direction = 1) +
  theme_minimal() +
  labs(x = "Number of Phenotypes",
       y = "Mean of Computation Time (Second)",
       title="Mean of Computational Time for Different Number of Subjects and Phenotypes",
       subtitle = "1,000 simulations are applied") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10)) +
  theme(plot.title = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) 
dev.off()



############################# 500 Simulation #############################
rm(list = ls())
## ---- 500 subjects + 500 simulations ---- ###
all_files <- list.files("Result/Oct52021GaMuT468Check/")
subject_500_dglm_2000 <- grep("*._GAMuT_time_adj_500_subject_dglm_500.*", all_files, value = T)

for (i in 1:length(subject_500_dglm_2000)) {
  file_name <- subject_500_dglm_2000[i]
  load(paste("Result/Oct52021GaMuT468Check/", file_name, sep = ""))
  pheno = substr(file_name, nchar(file_name) - 7, nchar(file_name) - 6)
  assign(paste("subject_500_dglm_2000", pheno, "_Pheno", sep = ""), GAMuT_time)
}

CCT_pvalue_df_1 <- data.frame(Number_Sub = 500,
                              Number_Sim = 500,
                              subject_500_dglm_2000_4_Pheno,
                              subject_500_dglm_2000_6_Pheno,
                              subject_500_dglm_2000_8_Pheno)
apply(CCT_pvalue_df_1[, 3:5], 2, mean)
# subject_500_dglm_2000_4_Pheno subject_500_dglm_2000_6_Pheno subject_500_dglm_2000_8_Pheno 
#                    0.4812881                     0.5098423                     0.5295877 
names(CCT_pvalue_df_1)[3:5] <- c("4_pheno", "6_pheno", "8_pheno")

## ---- 1000 subjects + 2000 simulations ---- ###
all_files <- list.files("Result/Oct52021GaMuT468Check/")
subject_1000_dglm_2000 <- grep("*._GAMuT_time_adj_1000_subject_dglm_500.*", all_files, value = T)

for (i in 1:length(subject_1000_dglm_2000)) {
  file_name <- subject_1000_dglm_2000[i]
  load(paste("Result/Oct52021GaMuT468Check/", file_name, sep = ""))
  pheno = substr(file_name, nchar(file_name) - 7, nchar(file_name) - 6)
  assign(paste("subject_1000_dglm_2000", pheno, "_Pheno", sep = ""), GAMuT_time)
}

CCT_pvalue_df_2 <- data.frame(Number_Sub = 1000,
                              Number_Sim = 500,
                              subject_1000_dglm_2000_4_Pheno,
                              subject_1000_dglm_2000_6_Pheno,
                              subject_1000_dglm_2000_8_Pheno)
apply(CCT_pvalue_df_2[, 3:5], 2, mean)
# subject_1000_dglm_2000_4_Pheno subject_1000_dglm_2000_6_Pheno subject_1000_dglm_2000_8_Pheno 
#                       3.444421                       3.445385                       3.482219 
names(CCT_pvalue_df_2)[3:5] <- c("4_pheno", "6_pheno", "8_pheno")

df_2000_simu <- rbind(CCT_pvalue_df_1, 
                      CCT_pvalue_df_2) 

temp = 
  df_2000_simu %>%
  dplyr::select(., -Number_Sim) %>%
  gather(., "Number Phenotypes", "Computation Time", -Number_Sub)%>%
  mutate(Number_SubF = factor(Number_Sub))
table(temp$Number_Sub, temp$Number_SubF)

# ---- Median
ggsave("Analysis//Oct52021GaMuT468Check/median_computationalTime_over_numberPheno_500simu.pdf", width = 12, height = 6)
df_2000_simu %>%
  dplyr::select(., -Number_Sim) %>%
  gather(., "Number Phenotypes", "Computation Time", - Number_Sub)%>%
  mutate("Number Subjects" = factor(Number_Sub)) %>%
  dplyr::select(., -Number_Sub) %>%
  group_by(`Number Subjects`, `Number Phenotypes`) %>%
  summarise(Median = median(`Computation Time`)) %>%
  ggplot(., aes(x=`Number Phenotypes`, y=Median, fill= `Number Subjects`)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label = round(Median,2)), size = 3,  position = position_dodge(0.9)) +
  scale_fill_brewer(direction = 1) +
  theme_minimal() +
  labs(x = "Number of Phenotypes",
       y = "Median of Computation Time (Second)",
       title="Median of Computational Time for Different Number of Subjects and Phenotypes",
       subtitle = "500 simulations are applied") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10)) +
  theme(plot.title = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) 
dev.off()

# ---- Mean
ggsave("Analysis//Oct52021GaMuT468Check/mean_computationalTime_over_numberPheno_500simu.pdf", width = 12, height = 6)
df_2000_simu %>%
  dplyr::select(., -Number_Sim) %>%
  gather(., "Number Phenotypes", "Computation Time", - Number_Sub)%>%
  mutate("Number Subjects" = factor(Number_Sub)) %>%
  dplyr::select(., -Number_Sub) %>%
  group_by(`Number Subjects`, `Number Phenotypes`) %>%
  summarise(Mean = mean(`Computation Time`)) %>%
  ggplot(., aes(x=`Number Phenotypes`, y=Mean, fill= `Number Subjects`)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label = round(Mean,2)), size = 3,  position = position_dodge(0.9)) +
  scale_fill_brewer(direction = 1) +
  theme_minimal() +
  labs(x = "Number of Phenotypes",
       y = "Mean of Computation  Time (Second)",
       title="Mean of Computational Time for Different Number of Subjects and Phenotypes",
       subtitle = "500 simulations are applied") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10)) +
  theme(plot.title = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) 
dev.off()


