#!/bin/bash

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=ALL
### please enter your own email address below in order to track the results
#SBATCH --mail-user=apatte12@uwyo.edu
### enter any job name that you prefer
#SBATCH --job-name=MI_1
#SBATCH --array=1-7500


module load arcc/1.0 gcc/14.2.0 r/4.4.0

cd /project/modelscape/analyses/MissingTS/missing-data

config=Model_Runs/poiss_Sim_mods/RickerConfigA_6.txt

datFile=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
parFile=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
clsize=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
saveFile=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)
index1=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)
index2=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $config)
seed=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $config)

Rscript Model_Runs/poiss_Sim_mods/modelruns_ricker_MI.R ${datFile} ${parFile} ${clsize} ${saveFile} ${index1} ${index2} ${seed} > Model_Runs/poiss_trial/outputRickerTrial_${SLURM_ARRAY_TASK_ID}.txt




max_runs=20
cur_runs=1
while [ $cur_runs -lt $max_runs ]
do

if [ -e Model_Runs/poiss_Sim_mods/RickerA_resultTable_$SLURM_ARRAY_TASK_ID.csv ]
then
    echo "yay this run is done"
    cur_runs=100
else
    echo "shoot trying to rerun and current run number is: " $cur_runs
    newseed=$((${seed}+$cur_runs))
    Rscript Model_Runs/poiss_Sim_mods/modelruns_ricker_MI.R ${datFile} ${parFile} ${clsize} ${saveFile} ${index1} ${index2} $newseed $cur_runs > Model_Runs/poiss_trial/outputRickerTrial_${SLURM_ARRAY_TASK_ID}.txt
    ((cur_runs++))
fi

done
echo "Out of the loop"

