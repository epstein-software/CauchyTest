library(tidyverse)

setwd("~/Dropbox/Emory Courses/DOHG/")


###### Predefined global parameter
#' Number of simulations, number of subjects, Number of phenotypes and MAF
#' (independent of covariates)

Nsub_array   <- c(5000, 10000)
Npheno_array <- c(4, 6, 8, 10, 12)
maf_array    <- 0.25
gamma_array  <- c(NA, 0.15, 0.55, 0.95)
delta_array  <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.25, 0.35, 0.45, 0.6, 1.5) # add some more gamma
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
  
  file_name_SMAT_pvalue               <- paste("pvalue_power", "SMAT", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_name_includeXsquare_SMAT_pvalue <- paste("pvalue_power", "includeXsquare_SMAT", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_name_CCT_pvalue                <- paste("pvalue_power", "CCT", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")
  file_name_includeXsquare_CCT_pvalue <- paste("pvalue_power", "includeXsquare_CCT", "adj", "dglm", Nsim, "subj", Nsub, "pheno", Npheno, "MAF", maf, "Gamma", gamma, "Delta", delta, "DeltaSelect", delta_select, sep = "_")

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

gamma_array  <- c(NA, 0.15, 0.55, 0.95)
delta_array  <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.25, 0.35, 0.45, 0.6, 1.5) # add some more gamma

param_comb_plot1 <- expand.grid(gamma  = gamma_array,
                                delta  = delta_array)

for (index in 1:nrow(param_comb_plot1)) {
  gamma_val <- param_comb_plot1[index, "gamma"]
  delta_val <- param_comb_plot1[index, "delta"]
  
  if (is.na(gamma_val)){
    temp_df <- full_pvalue_df %>%
      select(-maf) %>%
      filter(is.na(gamma),
             delta == delta_val)    
  } else {
    temp_df <- full_pvalue_df %>%
      select(-maf) %>%
      filter(gamma == gamma_val,
             delta %in% delta_val)
  }
  
  if (nrow(temp_df) != 400000) {
    stop("data slice is incorrect")
  }
  
  temp_df_F <- temp_df %>%
    select(-gamma, -delta) %>%
    gather(., 
           key = "Method", 
           value = "Power",
           SMAT,
           SMAT_includeXsquare,
           CCT,
           CCT_includeXsquare) %>%
    mutate(less_than_alpha = ifelse(Power <= 0.05, 1, 0),
           less_than_alpha_one_th = ifelse(Power <= 0.05/10000, 1, 0),
           less_than_alpha_fiv_th = ifelse(Power <= 0.05/50000, 1, 0),
           less_than_alpha_ten_th = ifelse(Power <= 0.05/100000, 1, 0),
           less_than_alpha_eig_th = ifelse(Power <= 0.05/1000000, 1, 0),
           is_zero = ifelse(Power == 0, 1, 0)) %>%
    group_by(., 
             Nsub,
             Npheno,
             delta_select,
             Method) %>%
    summarise(Total_Count = n(),
              Percent_less_than_alpha = sum(less_than_alpha)/Total_Count,
              Percent_less_than_alpha_one_th = sum(less_than_alpha_one_th)/Total_Count,
              Percent_less_than_alpha_fiv_th = sum(less_than_alpha_fiv_th)/Total_Count,
              Percent_less_than_alpha_ten_th = sum(less_than_alpha_ten_th)/Total_Count,
              Percent_less_than_alpha_eig_th = sum(less_than_alpha_eig_th)/Total_Count,
              Percent_is_zero                = sum(is_zero)/Total_Count) 
  temp_df_F$delta_select_F = factor(temp_df_F$delta_select)
  
  # Un-comment this if you remove the one with 0 and produce the plots starting with
  # remove_zero_
  # temp_df_F = temp_df_F[temp_df_F$Percent_is_zero!=1,]
  if (all(temp_df_F$Percent_less_than_alpha - temp_df_F$Percent_less_than_alpha_one_th>= 0) == FALSE){
    stop ("error with the alpha and alpha 2")
  }
  
  if (all(temp_df_F$Percent_less_than_alpha_one_th - temp_df_F$Percent_less_than_alpha_fiv_th>= 0) == FALSE){
    stop ("error with the alpha 2 and alpha 3")
  }
  
  if (all(temp_df_F$Percent_less_than_alpha_fiv_th - temp_df_F$Percent_less_than_alpha_ten_th>= 0) == FALSE){
    stop ("error with the alpha 3 and alpha 4")
  }
  
  ggplot(temp_df_F, 
         aes(x = delta_select_F, y = Percent_less_than_alpha,  group=Method)) +
    geom_line(aes(linetype=Method, color = Method))+
    ggtitle(paste("Power Simulation with gamma =", gamma_val, "and delta =", delta_val)) +
    labs(subtitle = expression(The~pre-specified~threshold~is~5~X~10^-2)) +
    xlab("Sparsity: decreases from 0 to 1") + 
    ylab(expression(Percent~less~than~the~threshold~5~X~10^-2)) +
    theme(
      plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
      axis.title.x = element_text(color="blue", size=16, face="bold"),
      axis.title.y = element_text(color="#993333", size=16, face="bold")
    ) +
    theme_bw() +
    facet_wrap(~ Nsub + Npheno, ncol = 5) + 
    theme(strip.text = element_text(size=12, face="bold", color="darkblue"))
  ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_delta_01242022/", index, "_thre_5_minus_2_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_val, ".pdf", sep = ""),
         width = 40, height = 20, units = "cm")
  
  # ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_delta/remove_zero_", index, "_thre_5_minus_2_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_val, ".pdf", sep = ""),
  #       width = 40, height = 20, units = "cm")
  
  
  ggplot(temp_df_F, 
         aes(x = delta_select_F, y = Percent_less_than_alpha_one_th,  group=Method)) +
    geom_line(aes(linetype=Method, color = Method))+
    ggtitle(paste("Power Simulation with gamma =", gamma_val, "and delta =", delta_val)) +
    labs(subtitle = expression(The~pre-specified~threshold~is~5~X~10^-6)) +
    xlab("Sparsity: decreases from 0 to 1") + 
    #ylab(expression("Percent less than the specified threshold: 0.05/10000", expression(5 X 10^{-6}))) +
    ylab(expression(Percent~less~than~the~threshold~5~X~10^-6)) +
    theme(
      plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
      axis.title.x = element_text(color="blue", size=16, face="bold"),
      axis.title.y = element_text(color="#993333", size=16, face="bold")
    ) +
    theme_bw() +
    facet_wrap(~ Nsub + Npheno, ncol = 5) + 
    theme(strip.text = element_text(size=12, face="bold", color="darkblue"))
  ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_delta_01242022/", index, "_thre_5_minus_6_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_val, ".pdf", sep = ""),
         width = 40, height = 20, units = "cm")
  # ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_delta/remove_zero_", index, "_thre_5_minus_6_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_val, ".pdf", sep = ""),
  #        width = 40, height = 20, units = "cm")
  
  
  ggplot(temp_df_F, 
         aes(x = delta_select_F, y = Percent_less_than_alpha_fiv_th,  group=Method)) +
    geom_line(aes(linetype=Method, color = Method))+
    ggtitle(paste("Power Simulation with gamma =", gamma_val, "and delta =", delta_val)) +
    labs(subtitle = expression(The~pre-specified~threshold~is~1~X~10^-6)) +
    xlab("Sparsity: decreases from 0 to 1") + 
    #ylab(expression("Percent less than the specified threshold: 0.05/10000", expression(5 X 10^{-6}))) +
    ylab(expression(Percent~less~than~the~threshold~1~X~10^-6)) +
    theme(
      plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
      axis.title.x = element_text(color="blue", size=16, face="bold"),
      axis.title.y = element_text(color="#993333", size=16, face="bold")
    ) +
    theme_bw() +
    facet_wrap(~ Nsub + Npheno, ncol = 5) + 
    theme(strip.text = element_text(size=12, face="bold", color="darkblue"))
  ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_delta_01242022/", index, "_thre_1_minus_6_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_val, ".pdf", sep = ""),
         width = 40, height = 20, units = "cm")
  # ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_delta/remove_zero_", index, "_thre_1_minus_6_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_val, ".pdf", sep = ""),
  #        width = 40, height = 20, units = "cm")
  
  
  ggplot(temp_df_F, 
         aes(x = delta_select_F, y = Percent_less_than_alpha_ten_th,  group=Method)) +
    geom_line(aes(linetype=Method, color = Method))+
    ggtitle(paste("Power Simulation with gamma =", gamma_val, "and delta =", delta_val)) +
    labs(subtitle = expression(The~pre-specified~threshold~is~1~X~10^-7)) +
    xlab("Sparsity: decreases from 0 to 1") + 
    #ylab(expression("Percent less than the specified threshold: 0.05/10000", expression(5 X 10^{-6}))) +
    ylab(expression(Percent~less~than~the~threshold~1~X~10^-7)) +
    theme(
      plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
      axis.title.x = element_text(color="blue", size=16, face="bold"),
      axis.title.y = element_text(color="#993333", size=16, face="bold")
    ) +
    theme_bw() +
    facet_wrap(~ Nsub + Npheno, ncol = 5) + 
    theme(strip.text = element_text(size=12, face="bold", color="darkblue"))
  ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_delta_01242022/", index, "_thre_1_minus_7_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_val, ".pdf", sep = ""),
         width = 40, height = 20, units = "cm")
  # ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_delta/remove_zero_", index, "_thre_1_minus_7_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_val, ".pdf", sep = ""),
  #        width = 40, height = 20, units = "cm")

  ggplot(temp_df_F, 
         aes(x = delta_select_F, y = Percent_less_than_alpha_eig_th,  group=Method)) +
    geom_line(aes(linetype=Method, color = Method))+
    ggtitle(paste("Power Simulation with gamma =", gamma_val, "and delta =", delta_val)) +
    labs(subtitle = expression(The~pre-specified~threshold~is~5~X~10^-8)) +
    xlab("Sparsity: decreases from 0 to 1") + 
    #ylab(expression("Percent less than the specified threshold: 0.05/1000000", expression(5 X 10^{-8}))) +
    ylab(expression(Percent~less~than~the~threshold~5~X~10^-8)) +
    theme(
      plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
      axis.title.x = element_text(color="blue", size=16, face="bold"),
      axis.title.y = element_text(color="#993333", size=16, face="bold")
    ) +
    theme_bw() +
    facet_wrap(~ Nsub + Npheno, ncol = 5) + 
    theme(strip.text = element_text(size=12, face="bold", color="darkblue"))
  ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_delta_01242022/", index, "_thre_5_minus_8_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_val, ".pdf", sep = ""),
         width = 40, height = 20, units = "cm")

}

############### Fix number of subject, gamma and sparsity #########################

gamma_array  <- c(NA, 0.15, 0.55, 0.95)
delta_select_array  <- c(0.25, 0.5, 0.75, 1)

param_comb_plot2 <- expand.grid(gamma         = gamma_array,
                                delta_select  = delta_select_array)

for (index in 1:nrow(param_comb_plot2)) {
  gamma_val <- param_comb_plot2[index, "gamma"]
  delta_select_val <- param_comb_plot2[index, "delta_select"]
  
  if (is.na(gamma_val)){
    temp_df <- full_pvalue_df %>%
      dplyr::select(-maf) %>%
      filter(is.na(gamma),
             delta_select == delta_select_val)    
  } else {
    temp_df <- full_pvalue_df %>%
      dplyr::select(-maf) %>%
      filter(gamma == gamma_val,
             delta_select == delta_select_val)
  }
  
  if (nrow(temp_df) != 1200000) {
    stop("data slice is incorrect")
  }
  
  temp_df_F <- temp_df %>%
    dplyr::select(-gamma, -delta_select) %>%
    gather(., 
           key = "Method", 
           value = "Power",
           SMAT,
           SMAT_includeXsquare,
           CCT,
           CCT_includeXsquare) %>%
    mutate(less_than_alpha = ifelse(Power <= 0.05, 1, 0),
           less_than_alpha_one_th = ifelse(Power <= 0.05/10000, 1, 0),
           less_than_alpha_fiv_th = ifelse(Power <= 0.05/50000, 1, 0),
           less_than_alpha_ten_th = ifelse(Power <= 0.05/100000, 1, 0),
           less_than_alpha_eig_th = ifelse(Power <= 0.05/1000000, 1, 0),
           is_zero = ifelse(Power == 0, 1, 0)) %>%
    group_by(., 
             Nsub,
             Npheno,
             delta,
             Method) %>%
    summarise(Total_Count = n(),
              Percent_less_than_alpha = sum(less_than_alpha)/Total_Count,
              Percent_less_than_alpha_one_th = sum(less_than_alpha_one_th)/Total_Count,
              Percent_less_than_alpha_fiv_th = sum(less_than_alpha_fiv_th)/Total_Count,
              Percent_less_than_alpha_ten_th = sum(less_than_alpha_ten_th)/Total_Count,
              Percent_less_than_alpha_eig_th = sum(less_than_alpha_eig_th)/Total_Count,
              Percent_is_zero                = sum(is_zero)/Total_Count) 
  # temp_df_F$delta_F = factor(temp_df_F$delta)
  
  # Un-comment this if you remove the one with 0 and produce the plots starting with
  # remove_zero_
  # temp_df_F = temp_df_F[temp_df_F$Percent_is_zero!=1,]
  if (all(temp_df_F$Percent_less_than_alpha - temp_df_F$Percent_less_than_alpha_one_th>= 0) == FALSE){
    stop ("error with the alpha and alpha 2")
  }
  
  if (all(temp_df_F$Percent_less_than_alpha_one_th - temp_df_F$Percent_less_than_alpha_fiv_th>= 0) == FALSE){
    stop ("error with the alpha 2 and alpha 3")
  }
  
  if (all(temp_df_F$Percent_less_than_alpha_fiv_th - temp_df_F$Percent_less_than_alpha_ten_th>= 0) == FALSE){
    stop ("error with the alpha 3 and alpha 4")
  }
  
  ggplot(temp_df_F, 
         aes(x = delta, y = Percent_less_than_alpha,  group=Method)) +
    geom_line(aes(linetype=Method, color = Method))+
    ggtitle(paste("Power Simulation with gamma =", gamma_val, "and sparsity =", delta_select_val)) +
    labs(subtitle = expression(The~pre-specified~threshold~is~5~X~10^-2)) +
    xlab("Delta: 0.1, 0.2, 0.3, 0.35, 0.4, 0.45 0.5, 0.6, 0.75, 1, 1.25, 1.5") + 
    ylab(expression(Percent~less~than~the~threshold~5~X~10^-2)) +
    theme(
      plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
      axis.title.x = element_text(color="blue", size=16, face="bold"),
      axis.title.y = element_text(color="#993333", size=16, face="bold")
    ) +
    theme_bw() +
    facet_wrap(~ Nsub + Npheno, ncol = 5) + 
    theme(strip.text = element_text(size=12, face="bold", color="darkblue"))
  ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_sparsity_01242022/", index, "_thre_5_minus_2_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_select_val, ".pdf", sep = ""),
         width = 40, height = 20, units = "cm")
  # ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_delta/remove_zero_", index, "_thre_5_minus_2_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_select_val, ".pdf", sep = ""),
  #       width = 40, height = 20, units = "cm")
  
  ggplot(temp_df_F, 
         aes(x = delta, y = Percent_less_than_alpha_one_th,  group=Method)) +
    geom_line(aes(linetype=Method, color = Method))+
    ggtitle(paste("Power Simulation with gamma =", gamma_val, "and sparsity =", delta_select_val)) +
    labs(subtitle = expression(The~pre-specified~threshold~is~5~X~10^-6)) +
    xlab("Delta: 0.1, 0.2, 0.3, 0.35, 0.4, 0.45 0.5, 0.6, 0.75, 1, 1.25, 1.5") + 
    #ylab(expression("Percent less than the specified threshold: 0.05/10000", expression(5 X 10^{-6}))) +
    ylab(expression(Percent~less~than~the~threshold~5~X~10^-6)) +
    theme(
      plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
      axis.title.x = element_text(color="blue", size=16, face="bold"),
      axis.title.y = element_text(color="#993333", size=16, face="bold")
    ) +
    theme_bw() +
    facet_wrap(~ Nsub + Npheno, ncol = 5) + 
    theme(strip.text = element_text(size=12, face="bold", color="darkblue"))
  ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_sparsity_01242022/", index, "_thre_5_minus_6_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_select_val, ".pdf", sep = ""),
         width = 40, height = 20, units = "cm")
  # ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_delta/remove_zero_", index, "_thre_5_minus_6_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_select_val, ".pdf", sep = ""),
  #        width = 40, height = 20, units = "cm")
  
  if (delta_select_val > 0.5) {
    ggplot(temp_df_F, 
           aes(x = delta, y = Percent_less_than_alpha_fiv_th,  group=Method)) +
      geom_line(aes(linetype=Method, color = Method))+
      ggtitle(paste("Power Simulation with gamma =", gamma_val, "and sparsity =", delta_select_val)) +
      labs(subtitle = expression(The~pre-specified~threshold~is~1~X~10^-6)) +
      xlab("Delta: 0.1, 0.2, 0.3, 0.35, 0.4, 0.45 0.5, 0.6, 0.75, 1, 1.25, 1.5") + 
      #ylab(expression("Percent less than the specified threshold: 0.05/10000", expression(5 X 10^{-6}))) +
      ylab(expression(Percent~less~than~the~threshold~1~X~10^-6)) +
      theme(
        plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
        axis.title.x = element_text(color="blue", size=16, face="bold"),
        axis.title.y = element_text(color="#993333", size=16, face="bold")
      ) +
      theme_bw() +
      xlim(0, 0.8) +
      facet_wrap(~ Nsub + Npheno, ncol = 5) + 
      theme(strip.text = element_text(size=12, face="bold", color="darkblue"))
    ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_sparsity_01242022/", index, "_thre_1_minus_6_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_select_val, ".pdf", sep = ""),
           width = 40, height = 20, units = "cm")
  } else {
    ggplot(temp_df_F, 
           aes(x = delta, y = Percent_less_than_alpha_fiv_th,  group=Method)) +
      geom_line(aes(linetype=Method, color = Method))+
      ggtitle(paste("Power Simulation with gamma =", gamma_val, "and sparsity =", delta_select_val)) +
      labs(subtitle = expression(The~pre-specified~threshold~is~1~X~10^-6)) +
      xlab("Delta: 0.1, 0.2, 0.3, 0.35, 0.4, 0.45 0.5, 0.6, 0.75, 1, 1.25, 1.5") + 
      #ylab(expression("Percent less than the specified threshold: 0.05/10000", expression(5 X 10^{-6}))) +
      ylab(expression(Percent~less~than~the~threshold~1~X~10^-6)) +
      theme(
        plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
        axis.title.x = element_text(color="blue", size=16, face="bold"),
        axis.title.y = element_text(color="#993333", size=16, face="bold")
      ) +
      theme_bw() +
      facet_wrap(~ Nsub + Npheno, ncol = 5) + 
      theme(strip.text = element_text(size=12, face="bold", color="darkblue"))
    ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_sparsity_01242022/", index, "_thre_1_minus_6_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_select_val, ".pdf", sep = ""),
           width = 40, height = 20, units = "cm")
  }
  # ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_sparsity_01182022/remove_zero_", index, "_thre_1_minus_6_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_select_val, ".pdf", sep = ""),
  #        width = 40, height = 20, units = "cm")
  
  ggplot(temp_df_F, 
         aes(x = delta, y = Percent_less_than_alpha_ten_th,  group=Method)) +
    geom_line(aes(linetype=Method, color = Method))+
    ggtitle(paste("Power Simulation with gamma =", gamma_val, "and sparsity =", delta_select_val)) +
    labs(subtitle = expression(The~pre-specified~threshold~is~1~X~10^-7)) +
    xlab("Delta: 0.1, 0.2, 0.3, 0.35, 0.4, 0.45 0.5, 0.6, 0.75, 1, 1.25, 1.5") + 
    #ylab(expression("Percent less than the specified threshold: 0.05/10000", expression(5 X 10^{-6}))) +
    ylab(expression(Percent~less~than~the~threshold~1~X~10^-7)) +
    theme(
      plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
      axis.title.x = element_text(color="blue", size=16, face="bold"),
      axis.title.y = element_text(color="#993333", size=16, face="bold")
    ) +
    theme_bw() +
    facet_wrap(~ Nsub + Npheno, ncol = 5) + 
    theme(strip.text = element_text(size=12, face="bold", color="darkblue"))
  ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_sparsity_01242022/", index, "_thre_1_minus_7_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_select_val, ".pdf", sep = ""),
         width = 40, height = 20, units = "cm")
  # ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_sparsity_01182022/remove_zero_", index, "_thre_1_minus_7_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_select_val, ".pdf", sep = ""),
  #        width = 40, height = 20, units = "cm")
  
  if (delta_select_val > 0.5) {
    ggplot(temp_df_F, 
           aes(x = delta, y = Percent_less_than_alpha_eig_th,  group=Method)) +
      geom_line(aes(linetype=Method, color = Method))+
      ggtitle(paste("Power Simulation with gamma =", gamma_val, "and sparsity =", delta_select_val)) +
      labs(subtitle = expression(The~pre-specified~threshold~is~5~X~10^-8)) +
      xlab("Delta: 0.1, 0.2, 0.3, 0.35, 0.4, 0.45 0.5, 0.6, 0.75, 1, 1.25, 1.5") + 
      #ylab(expression("Percent less than the specified threshold: 0.05/10000", expression(5 X 10^{-6}))) +
      ylab(expression(Percent~less~than~the~threshold~5~X~10^-8)) +
      theme(
        plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
        axis.title.x = element_text(color="blue", size=16, face="bold"),
        axis.title.y = element_text(color="#993333", size=16, face="bold")
      ) +
      theme_bw() +
      xlim(0, 0.8) +
      facet_wrap(~ Nsub + Npheno, ncol = 5) + 
      theme(strip.text = element_text(size=12, face="bold", color="darkblue"))
    ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_sparsity_01242022/", index, "_thre_5_minus_8_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_select_val, ".pdf", sep = ""),
           width = 40, height = 20, units = "cm")
  } else {
    ggplot(temp_df_F, 
           aes(x = delta, y = Percent_less_than_alpha_eig_th,  group=Method)) +
      geom_line(aes(linetype=Method, color = Method))+
      ggtitle(paste("Power Simulation with gamma =", gamma_val, "and sparsity =", delta_select_val)) +
      labs(subtitle = expression(The~pre-specified~threshold~is~5~X~10^-8)) +
      xlab("Delta: 0.1, 0.2, 0.3, 0.35, 0.4, 0.45 0.5, 0.6, 0.75, 1, 1.25, 1.5") + 
      #ylab(expression("Percent less than the specified threshold: 0.05/10000", expression(5 X 10^{-6}))) +
      ylab(expression(Percent~less~than~the~threshold~5~X~10^-8)) +
      theme(
        plot.subtitle = element_text(color="black", size=16, face="bold.italic"),
        axis.title.x = element_text(color="blue", size=16, face="bold"),
        axis.title.y = element_text(color="#993333", size=16, face="bold")
      ) +
      theme_bw() +
      facet_wrap(~ Nsub + Npheno, ncol = 5) + 
      theme(strip.text = element_text(size=12, face="bold", color="darkblue"))
    ggsave(paste("Analysis/Power_Simulation_Dec2021/plot/fix_gamma_sparsity_01242022/", index, "_thre_5_minus_8_", "power_simulation_gamma_", gamma_val, "_delta_delta_val_", delta_select_val, ".pdf", sep = ""),
           width = 40, height = 20, units = "cm")
  }

  
}