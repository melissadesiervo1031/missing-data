#!/bin/bash

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=ALL
### please enter your own email address below in order to track the results
#SBATCH --mail-user=apatte12@uwyo.edu
### enter any job name that you prefer
#SBATCH --job-name=rickerRerunA
#SBATCH --array=1-7500


module load arcc/1.0 gcc/12.2.0 r/4.4.0

cd /project/modelscape/analyses/MissingTS/missing-data

config=Model_Runs/poiss_Sim_mods/RickerConfigA.txt

datFile=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
parFile=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
clsize=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
saveFile=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)
index1=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)
index2=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $config)


Rscript Model_Runs/poiss_Sim_mods/modelruns_rickerALT.R ${datFile} ${parFile} ${clsize} ${saveFile} ${index1} ${index2} > Model_Runs/poiss_Sim_mods/outputRickerA_${oldName}.txt

