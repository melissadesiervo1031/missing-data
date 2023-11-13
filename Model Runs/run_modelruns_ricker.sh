#!/bin/bash

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
### please enter your own email address below in order to track the results
#SBATCH --mail-user=apatte12@uwyo.edu
### enter any job name that you prefer
#SBATCH --job-name=rickerTestRun


module load arcc/1.0 gcc/12.2.0 r/4.2.2

cd /project/modelscape/analyses/MissingTS/missing-data

Rscript "Model Runs/modelruns_ricker.R" "data/missingDatasets/pois_sim_randMiss_A.rds" "data/missingDatasets/pois_sim_params.rds" 3 "Model Runs/RickerA_resultTable.rds" > "Model Runs/outputRickerA.txt"