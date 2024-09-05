library(here)
library(stats)
library(Amelia)
library(R.utils)

# this is a function just to test the amelia error- a clone of the fit_ricker_MI without returning results
# it doesn't return results because it doesn't make sense to fit a ricker model to this kind of data
# so it will just return a success message if it works
source(here("ModelChecking/fit_gaussian_test.R"))

# read in the gaussian datasets
datasets=readRDS(here("data/missingDatasets/gauss_sim_randMiss_A.rds"))


# Run this to experience the "R session aborted" error from Amelia and have to restart your R session
# Method can be set to "both", "lags", "lead", and "neither"
for(i in 1:50){
  for(j in 1:15){
    y1=datasets[[i]]$y[[1]][1:60] # retrieve a single timeseries and trim to length 60
    fit_gaussian_test(y1,5,ameliatimeout=15, p2samelia=1,method="both") # try run the imputations
  }
}




