library(tidyverse)

setwd("~/Dropbox/Emory Courses/DOHG")
source("Analysis/result_analysis_source_code.R")


############################# QQ-PLOT Type I Error #############################

## --- CCT --- ## 

load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_CCT_adj_dglm_10000_subj_1000_pheno_4.RData")
pval_CCT_4 <- pval_CCT_adj_dglm

load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_CCT_adj_dglm_10000_subj_1000_pheno_6.RData")
pval_CCT_6 <- pval_CCT_adj_dglm

load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_CCT_adj_dglm_10000_subj_1000_pheno_8.RData")
pval_CCT_8 <- pval_CCT_adj_dglm

load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_CCT_adj_dglm_10000_subj_1000_pheno_10.RData")
pval_CCT_10 <- pval_CCT_adj_dglm

load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_CCT_adj_dglm_10000_subj_1000_pheno_12.RData")
pval_CCT_12 <- pval_CCT_adj_dglm

load("Result/c4_6_Simulation_HPC/2021-09-28_pvalue_CCT_adj_dglm_10000_subj_1000_pheno_14.RData")
pval_CCT_14 <- pval_CCT_adj_dglm

load("Result/c4_6_Simulation_HPC/2021-09-28_pvalue_CCT_adj_dglm_10000_subj_1000_pheno_16.RData")
pval_CCT_16 <- pval_CCT_adj_dglm

load("Result/c4_6_Simulation_HPC/2021-09-28_pvalue_CCT_adj_dglm_10000_subj_1000_pheno_18.RData")
pval_CCT_18 <- pval_CCT_adj_dglm

load("Result/c4_6_Simulation_HPC/2021-10-03_pvalue_CCT_adj_dglm_10000_subj_1000_pheno_20.RData")
pval_CCT_20 <- pval_CCT_adj_dglm

CCT_pvalue_df <- data.frame(pval_CCT_4,
                            pval_CCT_6,
                            pval_CCT_8,
                            pval_CCT_10,
                            pval_CCT_12,
                            pval_CCT_14,
                            pval_CCT_16,
                            pval_CCT_18,
                            pval_CCT_20)
# 4 phenotypes
pdf(file="Analysis/c4_6_Simulation_HPC/qqplot_CCT_pvalue_4pheno.pdf")
qqunif.plot(CCT_pvalue_df$pval_CCT_4, 
            main = "Type 1 Error of Cauchy Combination Test \n4 Phenotypes on 1,000 subjects \n10,000 Simulations With 95% Confidence Interval",
            col = "#012169")
dev.off()
pdf(file="Analysis/c4_6_Simulation_HPC/hist_CCT_pvalue_4pheno.pdf")
ggplot(CCT_pvalue_df, aes(x = pval_CCT_4)) +
  geom_histogram(bins = 30, color = "#012169", fill = "white") +
  labs(x = "p-value",
       y = "Count",
       title="Histogram of Type 1 Error from Cauchy Combination Test \n4 Phenotypes on 1,000 subjects with 10,000 Simulations", y="Value", x="Statistics") +
  theme_minimal() +
  theme(plot.title = element_text(face="bold")) 
dev.off()
prop.table(table(CCT_pvalue_df$pval_CCT_4 <= 0.05)) #0.0492 

# 6 phenotypes
pdf(file="Analysis/c4_6_Simulation_HPC/qqplot_CCT_pvalue_6pheno.pdf")
qqunif.plot(CCT_pvalue_df$pval_CCT_6, 
            main = "Type 1 Error of Cauchy Combination Test \n6 Phenotypes on 1,000 subjects \n10,000 Simulations With 95% Confidence Interval",
            col = "#012169")
dev.off()
pdf(file="Analysis/c4_6_Simulation_HPC/hist_CCT_pvalue_6pheno.pdf")
ggplot(CCT_pvalue_df, aes(x = pval_CCT_6)) +
  geom_histogram(bins = 30, color = "#012169", fill = "white") +
  labs(x = "p-value",
       y = "Count",
       title="Histogram of Type 1 Error from Cauchy Combination Test \n6 Phenotypes on 1,000 subjects with 10,000 Simulations", y="Value", x="Statistics") +
  theme_minimal() +
  theme(plot.title = element_text(face="bold")) 
dev.off()
prop.table(table(CCT_pvalue_df$pval_CCT_6 <= 0.05)) #0.0469 

# 8 phenotypes
pdf(file="Analysis/c4_6_Simulation_HPC/qqplot_CCT_pvalue_8pheno.pdf")
qqunif.plot(CCT_pvalue_df$pval_CCT_8, 
            main = "Type 1 Error of Cauchy Combination Test \n8 Phenotypes on 1,000 subjects \n10,000 Simulations With 95% Confidence Interval",
            col = "#012169")
dev.off()
pdf(file="Analysis/c4_6_Simulation_HPC/hist_CCT_pvalue_8pheno.pdf")
ggplot(CCT_pvalue_df, aes(x = pval_CCT_8)) +
  geom_histogram(bins = 30, color = "#012169", fill = "white") +
  labs(x = "p-value",
       y = "Count",
       title="Histogram of Type 1 Error from Cauchy Combination Test \n8 Phenotypes on 1,000 subjects with 10,000 Simulations", y="Value", x="Statistics") +
  theme_minimal() +
  theme(plot.title = element_text(face="bold")) 
dev.off()
prop.table(table(CCT_pvalue_df$pval_CCT_8 <= 0.05)) #0.052 

# 10 phenotypes
pdf(file="Analysis/c4_6_Simulation_HPC/qqplot_CCT_pvalue_10pheno.pdf")
qqunif.plot(CCT_pvalue_df$pval_CCT_10, 
            main = "Type 1 Error of Cauchy Combination Test \n10 Phenotypes on 1,000 subjects \n10,000 Simulations With 95% Confidence Interval",
            col = "#012169")
dev.off()
pdf(file="Analysis/c4_6_Simulation_HPC/hist_CCT_pvalue_10pheno.pdf")
ggplot(CCT_pvalue_df, aes(x = pval_CCT_10)) +
  geom_histogram(bins = 30, color = "#012169", fill = "white") +
  labs(x = "p-value",
       y = "Count",
       title="Histogram of Type 1 Error from Cauchy Combination Test \n10 Phenotypes on 1,000 subjects with 10,000 Simulations", y="Value", x="Statistics") +
  theme_minimal() +
  theme(plot.title = element_text(face="bold")) 
dev.off()
prop.table(table(CCT_pvalue_df$pval_CCT_10 <= 0.05)) #0.0479  

# 12 phenotypes
pdf(file="Analysis/c4_6_Simulation_HPC/qqplot_CCT_pvalue_12pheno.pdf")
qqunif.plot(CCT_pvalue_df$pval_CCT_12, 
            main = "Type 1 Error of Cauchy Combination Test \n12 Phenotypes on 1,000 subjects \n10,000 Simulations With 95% Confidence Interval",
            col = "#012169")
dev.off()
pdf(file="Analysis/c4_6_Simulation_HPC/hist_CCT_pvalue_12pheno.pdf")
ggplot(CCT_pvalue_df, aes(x = pval_CCT_10)) +
  geom_histogram(bins = 30, color = "#012169", fill = "white") +
  labs(x = "p-value",
       y = "Count",
       title="Histogram of Type 1 Error from Cauchy Combination Test \n12 Phenotypes on 1,000 subjects with 10,000 Simulations", y="Value", x="Statistics") +
  theme_minimal() +
  theme(plot.title = element_text(face="bold")) 
dev.off()
prop.table(table(CCT_pvalue_df$pval_CCT_12 <= 0.05)) #0.0454   

# 14 phenotypes
pdf(file="Analysis/c4_6_Simulation_HPC/qqplot_CCT_pvalue_14pheno.pdf")
qqunif.plot(CCT_pvalue_df$pval_CCT_14, 
            main = "Type 1 Error of Cauchy Combination Test \n14 Phenotypes on 1,000 subjects \n10,000 Simulations With 95% Confidence Interval",
            col = "#012169")
dev.off()
pdf(file="Analysis/c4_6_Simulation_HPC/hist_CCT_pvalue_14pheno.pdf")
ggplot(CCT_pvalue_df, aes(x = pval_CCT_14)) +
  geom_histogram(bins = 30, color = "#012169", fill = "white") +
  labs(x = "p-value",
       y = "Count",
       title="Histogram of Type 1 Error from Cauchy Combination Test \n14 Phenotypes on 1,000 subjects with 10,000 Simulations", y="Value", x="Statistics") +
  theme_minimal() +
  theme(plot.title = element_text(face="bold")) 
dev.off()
prop.table(table(CCT_pvalue_df$pval_CCT_14 <= 0.05)) #0.0494    

# 16 phenotypes
pdf(file="Analysis/c4_6_Simulation_HPC/qqplot_CCT_pvalue_16pheno.pdf")
qqunif.plot(CCT_pvalue_df$pval_CCT_16, 
            main = "Type 1 Error of Cauchy Combination Test \n16 Phenotypes on 1,000 subjects \n10,000 Simulations With 95% Confidence Interval",
            col = "#012169")
dev.off()
pdf(file="Analysis/c4_6_Simulation_HPC/hist_CCT_pvalue_16pheno.pdf")
ggplot(CCT_pvalue_df, aes(x = pval_CCT_16)) +
  geom_histogram(bins = 30, color = "#012169", fill = "white") +
  labs(x = "p-value",
       y = "Count",
       title="Histogram of Type 1 Error from Cauchy Combination Test \n16 Phenotypes on 1,000 subjects with 10,000 Simulations", y="Value", x="Statistics") +
  theme_minimal() +
  theme(plot.title = element_text(face="bold")) 
dev.off()
prop.table(table(CCT_pvalue_df$pval_CCT_16 <= 0.05)) #0.0471     


# 18 phenotypes
pdf(file="Analysis/c4_6_Simulation_HPC/qqplot_CCT_pvalue_18pheno.pdf")
qqunif.plot(CCT_pvalue_df$pval_CCT_18, 
            main = "Type 1 Error of Cauchy Combination Test \n18 Phenotypes on 1,000 subjects \n10,000 Simulations With 95% Confidence Interval",
            col = "#012169")
dev.off()
pdf(file="Analysis/c4_6_Simulation_HPC/hist_CCT_pvalue_18pheno.pdf")
ggplot(CCT_pvalue_df, aes(x = pval_CCT_18)) +
  geom_histogram(bins = 30, color = "#012169", fill = "white") +
  labs(x = "p-value",
       y = "Count",
       title="Histogram of Type 1 Error from Cauchy Combination Test \n18 Phenotypes on 1,000 subjects with 10,000 Simulations", y="Value", x="Statistics") +
  theme_minimal() +
  theme(plot.title = element_text(face="bold")) 
dev.off()
prop.table(table(CCT_pvalue_df$pval_CCT_18 <= 0.05)) #0.0462 


# 20 phenotypes
pdf(file="Analysis/c4_6_Simulation_HPC/qqplot_CCT_pvalue_20pheno.pdf")
qqunif.plot(CCT_pvalue_df$pval_CCT_20, 
            main = "Type 1 Error of Cauchy Combination Test \n20 Phenotypes on 1,000 subjects \n10,000 Simulations With 95% Confidence Interval",
            col = "#012169")
dev.off()
pdf(file="Analysis/c4_6_Simulation_HPC/hist_CCT_pvalue_20pheno.pdf")
ggplot(CCT_pvalue_df, aes(x = pval_CCT_20)) +
  geom_histogram(bins = 30, color = "#012169", fill = "white") +
  labs(x = "p-value",
       y = "Count",
       title="Histogram of Type 1 Error from Cauchy Combination Test \n20 Phenotypes on 1,000 subjects with 10,000 Simulations", y="Value", x="Statistics") +
  theme_minimal() +
  theme(plot.title = element_text(face="bold")) 
dev.off()
prop.table(table(CCT_pvalue_df$pval_CCT_20 <= 0.05)) #0.0453 

########### Computation time for a single function #############################
final_df <- data.frame()

# 4 phenotypes
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_CCT_time_adj_dglm_10000_subj_1000_pheno_4.RData")
CCT_time_4 <- CCT_time 
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_GAMuT_time_adj_dglm_10000_subj_1000_pheno_4.RData")
GAMuT_time_4 <- GAMuT_time 
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_SMAT_time_adj_dglm_10000_subj_1000_pheno_4.RData")
SMAT_time_4 <- SMAT_time 

stats_df = as.data.frame(cbind(CCT_time_4,
                               GAMuT_time_4,
                               SMAT_time_4))
stats_df_dev <- stats_df %>%
  tidyr::gather(., "Method", "Time")
stats_df_dev$Method <- gsub("_.*", "", stats_df_dev$Method)
stats_df_dev$Phenotype <- 4
final_df <- rbind(final_df, stats_df_dev)

# fun_mean <- function(x){
#   return(data.frame(y=median(x, na.rm = TRUE),label=paste0("Mean is ", round(mean(x, na.rm = TRUE), 2))))}
# pdf(file="Analysis/c4_6_Simulation_HPC/comp_time_3methods_pvalue_4pheno.pdf", width = 10, height = 7)
# ggplot(data=stats_df_dev, aes(x=Method , y=Time)) +
#   geom_boxplot() +
#   #scale_y_continuous(labels = scales::scientific) + 
#   stat_summary(fun.data = "fun_mean", geom="text", vjust=-1.5, hjust = 1.2, size=5) +
#   theme_bw() + 
#   labs(x = "Method",
#        y = "Computation Time (Second)",
#        title="Distribution of Computational Time for the Three Methods in Seconds \n1000 Simulations with 4 Phenotypes", y="Value", x="Statistics") +
#   theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10))
# dev.off()

# 6 phenotypes
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_CCT_time_adj_dglm_10000_subj_1000_pheno_6.RData")
CCT_time_6 <- CCT_time 
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_GAMuT_time_adj_dglm_10000_subj_1000_pheno_6.RData")
GAMuT_time_6 <- GAMuT_time 
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_SMAT_time_adj_dglm_10000_subj_1000_pheno_6.RData")
SMAT_time_6 <- SMAT_time 

stats_df = as.data.frame(cbind(CCT_time_6,
                               GAMuT_time_6,
                               SMAT_time_6))
stats_df_dev <- stats_df %>%
  tidyr::gather(., "Method", "Time")
stats_df_dev$Method <- gsub("_.*", "", stats_df_dev$Method)
stats_df_dev$Phenotype <- 6
final_df <- rbind(final_df, stats_df_dev)

# fun_mean <- function(x){
#   return(data.frame(y=median(x, na.rm = TRUE),label=paste0("Mean is ", round(mean(x, na.rm = TRUE), 2))))}
# pdf(file="Analysis/c4_6_Simulation_HPC/comp_time_3methods_pvalue_6pheno.pdf", width = 10, height = 7)
# ggplot(data=stats_df_dev, aes(x=Method , y=Time)) +
#   geom_boxplot() +
#   #scale_y_continuous(labels = scales::scientific) + 
#   stat_summary(fun.data = "fun_mean", geom="text", vjust=-1.5, hjust = 1.2, size=5) +
#   theme_bw() + 
#   labs(x = "Method",
#        y = "Computation Time (Second)",
#        title="Distribution of Computational Time for the Three Methods in Seconds \n1000 Simulations with 6 Phenotypes", y="Value", x="Statistics") +
#   theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10))
# dev.off()

# 8 phenotypes
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_CCT_time_adj_dglm_10000_subj_1000_pheno_8.RData")
CCT_time_8 <- CCT_time 
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_GAMuT_time_adj_dglm_10000_subj_1000_pheno_8.RData")
GAMuT_time_8 <- GAMuT_time 
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_SMAT_time_adj_dglm_10000_subj_1000_pheno_8.RData")
SMAT_time_8 <- SMAT_time 

stats_df = as.data.frame(cbind(CCT_time_8,
                               GAMuT_time_8,
                               SMAT_time_8))
stats_df_dev <- stats_df %>%
  tidyr::gather(., "Method", "Time")
stats_df_dev$Method <- gsub("_.*", "", stats_df_dev$Method)
stats_df_dev$Phenotype <- 8
final_df <- rbind(final_df, stats_df_dev)

# fun_mean <- function(x){
#   return(data.frame(y=median(x, na.rm = TRUE),label=paste0("Mean is ", round(mean(x, na.rm = TRUE), 2))))}
# pdf(file="Analysis/c4_6_Simulation_HPC/comp_time_3methods_pvalue_8pheno.pdf", width = 10, height = 7)
# ggplot(data=stats_df_dev, aes(x=Method , y=Time)) +
#   geom_boxplot() +
#   #scale_y_continuous(labels = scales::scientific) + 
#   stat_summary(fun.data = "fun_mean", geom="text", vjust=-1.5, hjust = 1.2, size=5) +
#   theme_bw() + 
#   labs(x = "Method",
#        y = "Computation Time (Second)",
#        title="Distribution of Computational Time for the Three Methods in Seconds \n1000 Simulations with 8 Phenotypes", y="Value", x="Statistics") +
#   theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10))
# dev.off()

# 10 phenotypes
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_CCT_time_adj_dglm_10000_subj_1000_pheno_10.RData")
CCT_time_10 <- CCT_time 
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_GAMuT_time_adj_dglm_10000_subj_1000_pheno_10.RData")
GAMuT_time_10 <- GAMuT_time 
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_SMAT_time_adj_dglm_10000_subj_1000_pheno_10.RData")
SMAT_time_10 <- SMAT_time 

stats_df = as.data.frame(cbind(CCT_time_10,
                               GAMuT_time_10,
                               SMAT_time_10))
stats_df_dev <- stats_df %>%
  tidyr::gather(., "Method", "Time")
stats_df_dev$Method <- gsub("_.*", "", stats_df_dev$Method)
stats_df_dev$Phenotype <- 10
final_df <- rbind(final_df, stats_df_dev)

# fun_mean <- function(x){
#   return(data.frame(y=median(x, na.rm = TRUE),label=paste0("Mean is ", round(mean(x, na.rm = TRUE), 2))))}
# pdf(file="Analysis/c4_6_Simulation_HPC/comp_time_3methods_pvalue_10pheno.pdf", width = 10, height = 7)
# ggplot(data=stats_df_dev, aes(x=Method , y=Time)) +
#   geom_boxplot() +
#   #scale_y_continuous(labels = scales::scientific) + 
#   stat_summary(fun.data = "fun_mean", geom="text", vjust=-1.5, hjust = 1.2, size=5) +
#   theme_bw() + 
#   labs(x = "Method",
#        y = "Computation Time (Second)",
#        title="Distribution of Computational Time for the Three Methods in Seconds \n1000 Simulations with 10 Phenotypes", y="Value", x="Statistics") +
#   theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10))
# dev.off()

# 12 phenotypes
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_CCT_time_adj_dglm_10000_subj_1000_pheno_12.RData")
CCT_time_12 <- CCT_time 
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_GAMuT_time_adj_dglm_10000_subj_1000_pheno_12.RData")
GAMuT_time_12 <- GAMuT_time 
load("Result/c4_6_Simulation_HPC/2021-09-27_pvalue_SMAT_time_adj_dglm_10000_subj_1000_pheno_12.RData")
SMAT_time_12 <- SMAT_time 

stats_df = as.data.frame(cbind(CCT_time_12,
                               GAMuT_time_12,
                               SMAT_time_12))
stats_df_dev <- stats_df %>%
  tidyr::gather(., "Method", "Time")
stats_df_dev$Method <- gsub("_.*", "", stats_df_dev$Method)
stats_df_dev$Phenotype <- 12
final_df <- rbind(final_df, stats_df_dev)

# fun_mean <- function(x){
#   return(data.frame(y=median(x, na.rm = TRUE),label=paste0("Mean is ", round(mean(x, na.rm = TRUE), 2))))}
# pdf(file="Analysis/c4_6_Simulation_HPC/comp_time_3methods_pvalue_12pheno.pdf", width = 10, height = 7)
# ggplot(data=stats_df_dev, aes(x=Method , y=Time)) +
#   geom_boxplot() +
#   #scale_y_continuous(labels = scales::scientific) + 
  # stat_summary(fun.data = "fun_mean", geom="text", vjust=-1.5, hjust = 1.2, size=5) +
  # theme_bw() +
  # labs(x = "Method",
  #      y = "Computation Time (Second)",
  #      title="Distribution of Computational Time for the Three Methods in Seconds \n1000 Simulations with 12 Phenotypes", y="Value", x="Statistics") +
  # theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10))
# dev.off()

# 14 phenotypes
load("Result/c4_6_Simulation_HPC/2021-09-28_pvalue_CCT_time_adj_dglm_10000_subj_1000_pheno_14.RData")
CCT_time_14 <- CCT_time 
load("Result/c4_6_Simulation_HPC/2021-09-28_pvalue_GAMuT_time_adj_dglm_10000_subj_1000_pheno_14.RData")
GAMuT_time_14 <- GAMuT_time 
load("Result/c4_6_Simulation_HPC/2021-09-28_pvalue_SMAT_time_adj_dglm_10000_subj_1000_pheno_14.RData")
SMAT_time_14 <- SMAT_time 

stats_df = as.data.frame(cbind(CCT_time_14,
                               GAMuT_time_14,
                               SMAT_time_14))
stats_df_dev <- stats_df %>%
  tidyr::gather(., "Method", "Time")
stats_df_dev$Method <- gsub("_.*", "", stats_df_dev$Method)
stats_df_dev$Phenotype <- 14
final_df <- rbind(final_df, stats_df_dev)

# 16 phenotypes
load("Result/c4_6_Simulation_HPC/2021-09-28_pvalue_CCT_time_adj_dglm_10000_subj_1000_pheno_16.RData")
CCT_time_16 <- CCT_time 
load("Result/c4_6_Simulation_HPC/2021-09-28_pvalue_GAMuT_time_adj_dglm_10000_subj_1000_pheno_16.RData")
GAMuT_time_16 <- GAMuT_time 
load("Result/c4_6_Simulation_HPC/2021-09-28_pvalue_SMAT_time_adj_dglm_10000_subj_1000_pheno_16.RData")
SMAT_time_16 <- SMAT_time 

stats_df = as.data.frame(cbind(CCT_time_16,
                               GAMuT_time_16,
                               SMAT_time_16))
stats_df_dev <- stats_df %>%
  tidyr::gather(., "Method", "Time")
stats_df_dev$Method <- gsub("_.*", "", stats_df_dev$Method)
stats_df_dev$Phenotype <- 16
final_df <- rbind(final_df, stats_df_dev)

# 18 phenotypes
load("Result/c4_6_Simulation_HPC/2021-09-28_pvalue_CCT_time_adj_dglm_10000_subj_1000_pheno_18.RData")
CCT_time_18 <- CCT_time 
load("Result/c4_6_Simulation_HPC/2021-09-28_pvalue_GAMuT_time_adj_dglm_10000_subj_1000_pheno_18.RData")
GAMuT_time_18 <- GAMuT_time 
load("Result/c4_6_Simulation_HPC/2021-09-28_pvalue_SMAT_time_adj_dglm_10000_subj_1000_pheno_18.RData")
SMAT_time_18 <- SMAT_time 

stats_df = as.data.frame(cbind(CCT_time_18,
                               GAMuT_time_18,
                               SMAT_time_18))
stats_df_dev <- stats_df %>%
  tidyr::gather(., "Method", "Time")
stats_df_dev$Method <- gsub("_.*", "", stats_df_dev$Method)
stats_df_dev$Phenotype <- 18
final_df <- rbind(final_df, stats_df_dev)

# 20 phenotypes
load("Result/c4_6_Simulation_HPC/2021-10-03_pvalue_CCT_time_adj_dglm_10000_subj_1000_pheno_20.RData")
CCT_time_20 <- CCT_time 
load("Result/c4_6_Simulation_HPC/2021-10-03_pvalue_GAMuT_time_adj_dglm_10000_subj_1000_pheno_20.RData")
GAMuT_time_20 <- GAMuT_time 
load("Result/c4_6_Simulation_HPC/2021-10-03_pvalue_SMAT_time_adj_dglm_10000_subj_1000_pheno_20.RData")
SMAT_time_20 <- SMAT_time 

stats_df = as.data.frame(cbind(CCT_time_20,
                               GAMuT_time_20,
                               SMAT_time_20))
stats_df_dev <- stats_df %>%
  tidyr::gather(., "Method", "Time")
stats_df_dev$Method <- gsub("_.*", "", stats_df_dev$Method)
stats_df_dev$Phenotype <- 20
final_df <- rbind(final_df, stats_df_dev)


# Plot
final_df$Phenotype <- factor(final_df$Phenotype)

ggsave("Analysis/c4_6_Simulation_HPC/computationalTime_over_numberPheno.pdf", width = 12, height = 6)
final_df %>%
  group_by(Phenotype, Method) %>%
  summarise(Median = median(Time)) %>%
  ggplot(., aes(x=Phenotype, y=Median, fill=Method)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(direction = -1) +
  theme_minimal() +
  labs(x = "Number of Phenotypes",
       y = "Median of Computation Time (Second)",
       title="Median of Computational Time for the Three Methods in Seconds",
       subtitle = "10,000 Simulations with 1,000 Subjects in Each Simulation for Different Number of Phenotypes") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10)) +
  theme(plot.title = element_text(size=14,face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) 
dev.off()
################ Overall Computation time for three function ###################
