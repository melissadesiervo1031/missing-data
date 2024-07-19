#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=26
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=Bias_checks_expanded.R
LogFile=bias_checks.log

# Change to the relevant working directory
cd /project/modelscape/analyses/MissingTS/missing-data/ModelChecking

# Load R
module load arcc/1.0  gcc/12.2.0 r/4.4.0

Rscript --vanilla $Rscript > $LogFile

