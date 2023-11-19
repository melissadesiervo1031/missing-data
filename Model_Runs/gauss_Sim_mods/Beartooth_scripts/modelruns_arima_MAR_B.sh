#!/bin/bash -l
#SBATCH --job-name modelruns_arima_MAR_B
#SBATCH --time=010:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=astears@uwyo.edu
#SBATCH --account=modelscape
#SBATCH --output=output/modelruns_arima_MAR_A_ex01_%B_%a.out
#SBATCH --array=0-4999

echo "SLURM_JOB_ID:" $SLURM_JOB_ID
echo "SLURM_JOB_NAME:" $SLURM_JOB_NAME
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "SLURM_ARRAY_TASK_ID:" $SLURM_ARRAY_TASK_ID

module load gcc/12.2.0 arcc/1.0 r/4.2.2

cd /project/modelscape/users/astears

Rscript --vanilla modelruns_arima_MAR_B.R $SLURM_ARRAY_TASK_ID

