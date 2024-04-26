#!/bin/bash

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --time=124:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=apatte12@uwyo.edu
#SBATCH --job-name=DA_stepsize
#SBATCH --array=1-20


module load arcc/1.0 gcc/12.2.0 r/4.2.2

cd /project/modelscape/analyses/MissingTS/missing-data

Rscript Model_Runs/DA_stepsize.R ${SLURM_ARRAY_TASK_ID} > Model_Runs/DA_stepsize_${SLURM_ARRAY_TASK_ID}.txt

