#!/bin/bash
#$ -cwd
#$ -N c4_6_sim_more_subj_larger_simulation
#$ -o /mnt/EpsteinFSS/data/sbian/CauchyTest/Result/New_4_6_Simulation/c4_6_sim_larger_simulation.$TASK_ID.out
#$ -e /mnt/EpsteinFSS/data/sbian/CauchyTest/Result/New_4_6_Simulation/c4_6_sim_larger_simulation.$TASK_ID.err
#$ -m e
#$ -t 1-75
#$ -tc 75

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
Rscript /mnt/EpsteinFSS/data/sbian/CauchyTest/Code/New_4_6_Simulation/c4_6_sim_more_subj.R -e "options(width = 120); sessioninfo::session_info()" --args ${SGE_TASK_ID}

echo "**** Job ends ****"
date