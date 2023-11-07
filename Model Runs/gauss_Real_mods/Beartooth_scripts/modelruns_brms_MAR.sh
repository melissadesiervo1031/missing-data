#!/bin/bash -l
#SBATCH --job-name modelruns_brms_gaussianRealMAR
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=astears@uwyo.edu
#SBATCH --account=modelscape
#SBATCH --output=modelruns_brms_MAR_gaussian_real_ex01_%A.out
#SBATCH --array=0-2
#SBATCH --mem=4G           # total memory per node

echo "SLURM_JOB_ID:" $SLURM_JOB_ID
echo "SLURM_JOB_NAME:" $SLURM_JOB_NAME
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "SLURM_ARRAY_TASK_ID:" $SLURM_ARRAY_TASK_ID

module load gcc/12.2.0 arcc/1.0 r/4.2.2

cd /project/modelscape/users/astears/

Rscript --vanilla modelruns_brms_MAR.R $SLURM_ARRAY_TASK_ID

