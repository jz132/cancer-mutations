#!/bin/bash
#
#SBATCH -J ana_promoter        # Job name
#SBATCH -n 1                     # Number of cores
#SBATCH --mem=16G             # Memory in MB
#SBATCH -o R_%A_%a.out
#SBATCH --mail-user=jz132@duke.edu  # Email where notifications will be sent
#SBATCH --mail-type=FAIL

#SBATCH --array=1-667%50

#Your actual work goes after this line

module load R/3.6.1-gcb02
Rscript analytical_promoter_fixed_new-slurm.R ${SLURM_ARRAY_TASK_ID} "BRCA-EU"
