#!/bin/bash

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
### please enter your own email address below in order to track the results
#SBATCH --mail-user=apatte12@uwyo.edu
### enter any job name that you prefer
#SBATCH --job-name=Brerun2


module load arcc/1.0 gcc/14.2.0 r/4.4.0

cd /project/modelscape/analyses/MissingTS/missing-data


Rscript Model_Runs/poiss_Sim_mods/modelruns_ricker.R data/missingDatasets/pois_sim_randMiss_B.rds data/missingDatasets/pois_sim_params.rds 5 Model_Runs/RickerB_randMiss_resultTable1.rds 1 2500 > Model_Runs/poiss_Sim_mods/outputRickerB_randMiss1.txt


