#!/bin/bash
#$ -cwd
#$ -N type_i_err_100000_sim
#$ -o /mnt/EpsteinFSS/data/sbian/CauchyTest/Result/New_4_6_Simulation_ExcludeGAMUT/New_4_6_Simulation_ExcludeGAMUT_type_i.$TASK_ID.out
#$ -e /mnt/EpsteinFSS/data/sbian/CauchyTest/Result/New_4_6_Simulation_ExcludeGAMUT/New_4_6_Simulation_ExcludeGAMUT_type_i.$TASK_ID.err
#$ -m e
#$ -t 1-40
#$ -tc 40

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
# module list

Rscript /mnt/EpsteinFSS/data/sbian/CauchyTest/Code/New_4_6_Simulation_ExcludeGAMUT/c4_6_sim.R -e "options(width = 120); sessioninfo::session_info()" --args ${SGE_TASK_ID}

echo "**** Job ends ****"
date