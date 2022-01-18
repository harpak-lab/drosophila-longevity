#!/bin/bash

#SBATCH -J int_model
#SBATCH -o OutputFiles/int_model.o%j
#SBATCH -e ErrorFiles/int_model.e%j
#SBATCH -p flex
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=rknowlton@utexas.edu
#SBATCH --array=1-272

module load Rstats
Rscript int_model.R $SLURM_ARRAY_TASK_ID

