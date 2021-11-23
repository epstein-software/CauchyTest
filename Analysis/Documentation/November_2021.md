## November 22, 2021

New Updates:

- Updated the slide *Presemtation/c46\\_Simulation\\_HPC/c4\\_6\\_v06.pptx*

- The computational time is analyzed using */Analysis/New_4_6_Power_Simulation/c4_6_sim_power_fix_delta_SMAT_CCT_LocalJob/c4_6_sim_power_fix_delta_LocalJob_OneSim_Analysis.R*
	+ The results were simulated using *c4_6_sim_power_random_select_fix_delta_LocalJob_OneSim.R*
	
```r
Nsub_array   <- c(500, 1000, 2000, 5000, 10000)
Npheno_array <- c(4, 6, 8, 10, 12)
maf_array    <- 0.25
gamma_array  <- c(0.15)
delta_array  <- c(0.5)
delta_select_array <- c(1)
Nsim         <- 4
```


- The power simulation is analyzed using *c4_6_sim_power_random_collect_fix_delta_SMAT_CCT_LocalJob_Analysis*
	+ The results were simulated using *c4_6_sim_power_random_select_fix_delta_LocalJob_OneSim.R*


```r
Nsub_array   <- c(500)
Npheno_array <- c(4, 6, 8, 10, 12)
maf_array    <- 0.25
gamma_array  <- c(0.15, 0.25, 0.35)
delta_array  <- c(0.5, 0.75, 1) 
delta_select_array <- c(0.25, 0.5, 0.75, 1)
Nsim         <- 1000
```

### Task this week:

- Look into the trick using eigen value for linear kernel used in GAMuT that might make the code to be more efficient.
- Reading the analysis result from real data from Dr. Wingo that has been transferred to HGCC
- Exploring the ROS/MAP transcriptomics gene expression.
- Exploring UKB data shared by Lori
- Getting some preliminary result for type I error. I might focus on SMAT and CCT. GAMuT needs a longer time to complete.
- Utilize GitHub
