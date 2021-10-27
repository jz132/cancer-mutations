#!/bin/bash

#SBATCH -o slurmlog/vm_%A_%a.out
#SBATCH -e slurmlog/vm_%A_%a.err

#SBATCH --array=1-667%67
#SBATCH --mem=2G
#SBATCH --mail-type=END
#SBATCH --mail-user=vm76@duke.edu

Rscript muteffect_entropy.R $SLURM_ARRAY_TASK_ID
