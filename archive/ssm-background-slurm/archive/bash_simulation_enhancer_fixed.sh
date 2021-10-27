#!/bin/bash
#
#SBATCH -J simulation        # Job name
#SBATCH -n 1                     # Number of cores
#SBATCH --mem=16G             # Memory in MB
#SBATCH -o R_%j.out              # File for STDOUT (with jobid = %j)
#SBATCH --mail-user=jz132@duke.edu  # Email where notifications will be sent
#SBATCH --mail-type=FAIL

#SBATCH --array=1-100

#Your actual work goes after this line

module load R/3.6.1-gcb02
Rscript simulation_enhancer_fixed-slurm-parallel.R ${SLURM_ARRAY_TASK_ID}
