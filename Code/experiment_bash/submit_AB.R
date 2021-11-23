#!/bin/bash
#$ -cwd
#$ -N null.dGLM.var
#$ -o logs/01-null/null.dGLM.$TASK_ID.txt
#$ -e logs/01-null/null.dGLM.$TASK_ID.txt
#$ -m e
#$ -t 1-9
#$ -tc 9

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
Rscript 01-null.R -e "options(width = 120); sessioninfo::session_info()" --args ${SGE_TASK_ID}

echo "**** Job ends ****"
date