###################################################################
# This script is designed to run on a cluster and takes
# arguments from the command line including:
#   1. The file storing the list of datasets.
#   2. The file storing the dataframe of the parameters
#     used to generate the data.
#   3. Cluster size for parallel computing.
#   4. file to which to save the results
#   5. Optional argument to specify the method to use (default is all).
# Results are stored in a tibble with simulation specs and a column
# with fitted values for each of the methods used to fit
# the models in the presence of missing data.
###################################################################

library(here)
library(parallel)
library(parallelly)
library(tidyverse)

# source functions
f_list <- list.files(here("Functions/"), full.names = T)
f_list=f_list[-grep("README",f_list)]
lapply(f_list, source)

# arguments from the shell
# Should include (1) data file (2) parameter file (3) number of nodes for cluster (4) save file location
# (5) beginnning index (6) ending index (7) optional- model list

# for testing outside of command line, can use this next line
#in_args=c("data/model_results/ricker_Sim_reruns/extinctSets.rds", "data/missingDatasets/pois_sim_params.rds", 2, "data/model_results/ricker_Sim_reruns/MI_reruns/extinctMI_1.csv", 1, 1,1493)
in_args <- commandArgs(trailingOnly = T)
cat(in_args)
# read in datafile
dat_flat <- readRDS(here(in_args[1]))
#pars <- readRDS(here(in_args[2]))
cat("data has been loaded")
set.seed(in_args[7])

# count number of missingness proportions
#nmiss_props <- length(dat[[1]]$y)

# get vector of actual proportion missing
prop_miss <- str_extract(
    names(dat_flat),
    pattern = "0.\\d+"
) %>% unlist() %>% as.numeric()


  
  # get input autocorrelation from names
  autocorrs <- str_extract(
    names(dat_flat),
    "autoCorr_\\d+"
  ) %>% str_extract(
    ., 
    pattern = "\\d+"
  ) %>% as.numeric() / 100
  
  # get actual autocorrelation from names
  autocorr_act <- str_extract(
      names(dat_flat),
      "Corr_0.\\d+"
    ) %>% str_extract(
      pattern = "0.\\d+"
    ) %>% as.numeric() %>% unlist()
  
  
  simNumber = as.integer(str_split(names(dat_flat), "_", simplify = TRUE)[,2] %>% 
                           str_extract("\\d+"))
  
  pars_part=cbind(simNumber,autocorrs,prop_miss,autocorr_act)
  


# time testing only
dat_flat <- dat_flat[as.numeric(in_args[5]):as.numeric(in_args[6])]
pars_part <- pars_part[as.numeric(in_args[5]):as.numeric(in_args[6]),]


### fitting the models ###

system.time({
  
  # store results in a vector
  results_vec <- numeric(8)
  names(results_vec) <- c("estim_r","estim_alpha","se_r","se_alpha","lower_r","lower_alpha","upper_r","upper_alpha")
  
  # only use MI
  
  res1 <- lapply(
    X = dat_flat,
    FUN = fit_ricker_MI_lead
  )
  
  if(is.na(res1[[1]][1][[1]][1])){
    results_vec[1]=res1[[1]][2][[1]][1] # estim_r
    results_vec[2]=NA # estim_alpha
    results_vec[3]=NA # se_r
    results_vec[4]=NA # se_alpha
    results_vec[5]=NA # lower_r
    results_vec[6]=NA # lower_alpha
    results_vec[7]=NA # upper_r
    results_vec[8]=NA # upper_alpha
  } else {
    # compile results
    results_vec[1]=res1[[1]][1][[1]][1] # estim_r
    results_vec[2]=res1[[1]][1][[1]][2] # estim_alpha
    results_vec[3]=res1[[1]][2][[1]][1] # se_r
    results_vec[4]=res1[[1]][2][[1]][2] # se_alpha
    results_vec[5]=res1[[1]][3][[1]][1] # lower_r
    results_vec[6]=res1[[1]][3][[1]][2] # lower_alpha
    results_vec[7]=res1[[1]][4][[1]][1] # upper_r
    results_vec[8]=res1[[1]][4][[1]][2] # upper_alpha
  }
  
})

if(is.na(in_args[8])){
  ran20=F
} else {
  ran20=as.numeric(in_args[8])>18
}
if(is.na(results_vec[2])&results_vec[1]!="population extinction"&!ran20){
  # don't save if result was NA, and reason wasn't extinction and we haven't run 20 yet
} else {
  # save results to file
  write.csv(cbind(t(pars_part),as.data.frame(t(results_vec))),file = here(in_args[4]))
}


