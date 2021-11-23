#!/bin/bash
#$ -cwd
#$ -N c4_6_sim_power_random_select_fix_delta_SMAT_CCT
#$ -o /mnt/EpsteinFSS/data/sbian/CauchyTest/Result/New_4_6_Power_Simulation/c4_6_sim_power_random_select_fix_delta_SMAT_CCT.$TASK_ID.out
#$ -e /mnt/EpsteinFSS/data/sbian/CauchyTest/Result/New_4_6_Power_Simulation/c4_6_sim_power_random_select_fix_delta_SMAT_CCT.$TASK_ID.err
#$ -m e
#$ -t 1-4500
#$ -tc 4500

echo "**** Job starts ****"
date

echo "**** HGCC info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module
module load R/4.1.1

## List current modules for reproducibility
module list

## Edit with your job command
Rscript /mnt/EpsteinFSS/data/sbian/CauchyTest/Code/New_4_6_Power_Simulation/c4_6_sim_power_random_select_fix_delta_SMAT_CCT.R -e "options(width = 120); sessioninfo::session_info()" --args ${SGE_TASK_ID}

echo "**** Job ends ****"
date