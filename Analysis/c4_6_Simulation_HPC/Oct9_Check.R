library(tidyverse)
setwd("~/Dropbox/Emory Courses/DOHG")
source("Analysis/result_analysis_source_code.R")

# This analysis is to double confirm the simulation for 4, 6, 8, 10 phenotypes

all_files <- list.files("Result/c4_6_Simulation_HPC/Oct9Check/")
GAMuT_Files <- grep("*.GAMuT_time_adj_dglm_10000_subj_1000_pheno.*", all_files, value = T)

for (i in 1:length(GAMuT_Files)) {
  file_name <- GAMuT_Files[i]
  load(paste("Result/c4_6_Simulation_HPC/Oct9Check/", file_name, sep = ""))
  pheno = substr(file_name, nchar(file_name) - 7, nchar(file_name) - 6)
  assign(paste("GAMuT_Files", pheno, "_Pheno", sep = ""), GAMuT_time)
}

CCT_pvalue_df_1 <- data.frame(Number_Sub = 1000,
                            Number_Sim = 10000,
                            GAMuT_Files_4_Pheno,
                            GAMuT_Files_6_Pheno,
                            GAMuT_Files_8_Pheno,
                            GAMuT_Files10_Pheno)
apply(CCT_pvalue_df_1[, 3:6], 2, mean)
apply(CCT_pvalue_df_1[, 3:6], 2, median)

names(CCT_pvalue_df_1)[3:6] <- c("4_pheno", "6_pheno", "8_pheno", "10_pheno")

df_2000_simu <- rbind(CCT_pvalue_df_1) 

temp = 
  df_2000_simu %>%
  dplyr::select(., -Number_Sim) %>%
  gather(., "Number Phenotypes", "Computation Time", -Number_Sub)%>%
  mutate(Number_SubF = factor(Number_Sub),
         `Number Phenotypes F` = factor(`Number Phenotypes`, levels = c("4_pheno", 
                                                                        "6_pheno",
                                                                        "8_pheno",
                                                                        "10_pheno")))
table(temp$Number_Sub, temp$Number_SubF)
table(temp$`Number Phenotypes F`, temp$`Number Phenotypes`)

# ---- Median
pdf("Analysis/c4_6_Simulation_HPC/recheck_median_computationalTime_first_four_pheno.pdf", width = 12, height = 6)
df_2000_simu %>%
  dplyr::select(., -Number_Sim) %>%
  gather(., "Number Phenotypes", "Computation Time", - Number_Sub)%>%
  mutate("Number Subjects" = factor(Number_Sub),
         `Number Phenotypes F` = factor(`Number Phenotypes`, levels = c("4_pheno", 
                                                                        "6_pheno",
                                                                        "8_pheno",
                                                                        "10_pheno"))) %>%
  dplyr::select(., -Number_Sub, -`Number Phenotypes`) %>%
  group_by(`Number Phenotypes F`) %>%
  summarise(Median = median(`Computation Time`)) %>%
  ggplot(., aes(x=`Number Phenotypes F`, y=Median)) + 
  geom_bar(stat="identity", position=position_dodge(), fill = "#1e376d") +
  geom_text(aes(label = round(Median,2)), size = 3,  vjust = -0.1) +
  theme_minimal() +
  labs(x = "Number of Phenotypes",
       y = "Median of Computation Time (Second)",
       title="Median of Computational Time for Different Number Phenotypes",
       subtitle = "10,000 simulations on 1,000 subjects are applied") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.3, size = 10)) +
  theme(plot.title = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) 
dev.off()

# ---- Mean
pdf("Analysis/c4_6_Simulation_HPC/recheck_mean_computationalTime_first_four_pheno.pdf", width = 12, height = 6)
df_2000_simu %>%
  dplyr::select(., -Number_Sim) %>%
  gather(., "Number Phenotypes", "Computation Time", - Number_Sub)%>%
  mutate("Number Subjects" = factor(Number_Sub),
         `Number Phenotypes F` = factor(`Number Phenotypes`, levels = c("4_pheno", 
                                                                        "6_pheno",
                                                                        "8_pheno",
                                                                        "10_pheno"))) %>%
  dplyr::select(., -Number_Sub, -`Number Phenotypes`) %>%
  group_by(`Number Phenotypes F`) %>%
  summarise(Mean = mean(`Computation Time`)) %>%
  ggplot(., aes(x=`Number Phenotypes F`, y=Mean)) + 
  geom_bar(stat="identity", fill = "#1e376d") +
  geom_text(aes(label = round(Mean,2)), size = 3,  vjust = -0.1) +
  theme_minimal() +
  labs(x = "Number of Phenotypes",
       y = "Mean of Computation Time (Second)",
       title="Mean of Computational Time for Different Number Phenotypes",
       subtitle = "10,000 simulations on 1,000 subjects are applied") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.3, size = 10)) +
  theme(plot.title = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) 
dev.off()


########### p-value
source("Analysis/result_analysis_source_code.R")

load("~/Dropbox/Emory Courses/DOHG/Result/c4_6_Simulation_HPC/Oct9Check/2021-10-09_pvalue_GAMuT_adj_dglm_10000_subj_1000_pheno_10.RData")
pdf("Analysis/c4_6_Simulation_HPC/GAMuT_pvalue_kernel.pdf", width = 12, height = 6)
hist(pval_GAMuT_adj_dglm[,1],
     main = "p-value for Phenotypic Similarity using Linear Kernel: \n10,000 simulations on 1,000 subjects \nare applied on 10 phenotypes")
prop.table(table(pval_GAMuT_adj_dglm[,1]<0.05)) # 0.0416
dev.off()

load("~/Dropbox/Emory Courses/DOHG/Result/c4_6_Simulation_HPC/Oct9Check/2021-10-09_pvalue_GAMuT_adj_dglm_10000_subj_1000_pheno_10.RData")
pdf("Analysis/c4_6_Simulation_HPC/GAMuT_pvalue_matrix.pdf", width = 12, height = 6)
hist(pval_GAMuT_adj_dglm[,2],
     main = "p-value fpr Phenotypic Similarity using Projection Matrix: \n10,000 simulations on 1,000 subjects \nare applied on 10 phenotypes")
prop.table(table(pval_GAMuT_adj_dglm[,2]<0.05)) # 0.042
dev.off()

# phenotypic similarity:  linear kernel
pdf("Analysis/c4_6_Simulation_HPC/GAMuT_qqplot_kernel.pdf", width = 12, height = 12)
qqunif.plot(pval_GAMuT_adj_dglm[,1], 
            main = "Type 1 Error of GAMuT Checking Phenotypic Similarity using Linear Kernel\n4 Phenotypes on 1,000 subjects \n10,000 Simulations With 95% Confidence Interval, 10 Phenotypes",
            col = "#012169")
dev.off()

# phenotypic similarity:  projection matrix
pdf("Analysis/c4_6_Simulation_HPC/GAMuT_qqplot_matrix.pdf", width = 12, height = 12)
qqunif.plot(pval_GAMuT_adj_dglm[,2], 
            main = "Type 1 Error of GAMuT Checking Phenotypic Similarity using Projection Matrix\n4 Phenotypes on 1,000 subjects \n10,000 Simulations With 95% Confidence Interval, 10 Phenotypes",
            col = "#012169")
dev.off()