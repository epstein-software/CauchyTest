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

## December 5th, 2021

This simulation will exclude GAMuT. Therefore, two new folders are created to perform this simulation

- New_4_6_Simulation_ExcludeGAMUT
	+ This is for type I error simulation


### Type I simulation

Your job-array 444755.1-40:1 ("type_i_err_100000_sim") has been submitted

#### Revision

1. Add X^2 to SMAT
2. Force the time to be seconds
3. Make the value of gamma a little bigger
4. Random seed is revised to ensure

The setting:

```r
Nsub_array   <- c(5000, 10000)
Npheno_array <- c(4, 6, 8, 10, 12)
maf_array    <- c(0.25)
gamma_array  <- c(NA, 0.15, 0.55, 0.95)
Nsim         <- 1
```

The seed might need to be modified as adding `set.seed(sim_Npheno*10 + array.id*2 + sample(1:500000000, 1))`.

```r

for (sim_Npheno in 1:Npheno) {
  set.seed(sim_Npheno*10 + array.id*2 + sample(1:500000000, 1))
  
  # randomly-generated intercept for phenotypes
  int_array[sim_Npheno] <- rnorm(1, 0, 5) 
  
  # randomly-generated main effect of SNP on phenotypes
  main_beta_g_array[sim_Npheno] <- runif(1, 0, 0.2)
  
  # randomly-generate main effect of continuous covariate on genotype
  # Random simulate gamma is gamma == NA
  if (is.na(gamma)) {
    main_gamma_con_array[sim_Npheno] <- runif(1, 0, 0.2)
  } else {
    main_gamma_con_array[sim_Npheno] <- gamma
  }
}
```

## December 19, 2021

### Type I simulation Update

For some unknown reason, HGCC deleted most of my jobs. There are currently only 5 jobs are running. Not sure what happened. 5 jobs have finished renning. I might need to re-submit the 30 jobs.

### Type II simulation

The script was under folder `Code\CNew_4_6_Power_Simulation_ExcludeGAMUT`. 

This job is running on Harvard cluster. HGCC does not have enough resource.

```r
Nsub_array   <- c(5000, 10000)
Npheno_array <- c(4, 6, 8, 10, 12)
maf_array    <- 0.25
gamma_array  <- c(NA, 0.15, 0.55, 0.95)
delta_array  <- c(0.5, 0.75, 1, 1.25)
delta_select_array <- c(0.25, 0.5, 0.75, 1)
Nsim         <- 10000



param_comb <- expand.grid(Nsub   = Nsub_array,
                          Npheno = Npheno_array,
                          maf    = maf_array,
                          gamma  = gamma_array,
                          delta  = delta_array,
                          delta_select = delta_select_array)
```

## December 19, 2021

### Type I simulation

```r
Nsub_array   <- c(5000, 10000)
Npheno_array <- c(4, 6, 8, 10, 12)
maf_array    <- c(0.25)
gamma_array  <- c(NA, 0.15, 0.55, 0.95)
job_slice_array <- seq(1, 100)
Nsim         <- 10000

```


## January 2, 2022

### Type I simulation

Type I simulation has completed. 

# January 9, 2022

### Type I simulation

Type I simulation has been completed. The plots have been plotted.

### Power simulation

Power simulation has completed. The plot has been completed.

Tried another set of delta to observe more detailed performance comparison between CCT and SMAT:

```r
delta_array  <- c(0.1, 0.2, 0.3, 0.4)

```


# January 14, 2022

Power simulation for 10 phenotypes only, c4_6_sim_power_20220113.R

Add more information: c4_6_sim_power_20220114.R

```r
Nsub_array   <- c(5000, 10000)
Npheno_array <- c(4, 6, 8, 10, 12)
maf_array    <- 0.25
gamma_array  <- c(NA, 0.15, 0.55, 0.95)
delta_array  <- c(0.35, 0.45, 0.6, 1.5)
delta_select_array <- c(0.25, 0.5, 0.75, 1)
Nsim         <- 10000

```


