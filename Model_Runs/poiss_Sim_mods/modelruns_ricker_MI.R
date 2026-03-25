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
#in_args=c("data/missingDatasets/pois_sim_randMiss_A.rds", "data/missingDatasets/pois_sim_params.rds", 2, "Model_Runs/RickerA_resultTable1.csv", 1648, 1648,1493)
#in_args=c("data/missingDatasets/pois_sim_randMiss_trim_A.rds", "data/missingDatasets/pois_sim_params.rds", 2, "Model_Runs/RickerA_resultTable1.csv", 142, 142,1493)
in_args <- commandArgs(trailingOnly = T)
cat(in_args)
# read in datafile
dat <- readRDS(here(in_args[1]))
pars <- readRDS(here(in_args[2]))
dat0=readRDS(here("data/ricker_0miss_datasets.rds"))
cat("data has been loaded")
set.seed(in_args[7])

# count number of missingness proportions
nmiss_props <- length(dat[[1]]$y)

# get vector of actual proportion missing
prop_miss <- map(
  dat,
  ~ str_extract(
    names(.x$y),
    pattern = "0.\\d+"
  )
) %>% unlist() %>% as.numeric()

# flatten the list
dat_flat <- map(
  dat,
  ~ pluck(.x, "y")
) %>% list_flatten()

dat_flat_o=dat_flat

if(!is.null(names(dat))){
  
  # get input autocorrelation from names
  autocorrs <- str_extract(
    names(dat),
    "autoCorr_\\d+"
  ) %>% str_extract(
    ., 
    pattern = "\\d+"
  ) %>% as.numeric() / 100
  
  # get actual autocorrelation from names
  autocorr_act <- map(
    dat,
    ~ str_extract(
      names(.x$y),
      "Corr_0.\\d+"
    ) %>% str_extract(
      pattern = "0.\\d+"
    ) %>% as.numeric()
  ) %>% unlist()
  
  # get number of different autocorrelations
  n_autocorrs <- length(unique(autocorrs))
  
  pars_full <- pars[rep(1:nrow(pars), n_autocorrs), ]
  pars_full <- pars_full[rep(1:nrow(pars_full), each = nmiss_props), ]
  
  pars_full <- pars_full %>% mutate(
    id = 1:length(dat_flat),
    autoCorr = rep(autocorrs, each = 16),
    propMiss = rep(
      seq(0.00, 0.75, by = 0.05),
      nrow(pars) * n_autocorrs
    ),
    actAutoCorr = autocorr_act,
    actPropMiss = prop_miss
  )
  
} else {
  
  pars_full <- pars[rep(1:nrow(pars), each = nmiss_props), ]
  pars_full <- pars_full %>% mutate(
    id = 1:length(dat_flat),
    propMiss = rep(
      seq(0.00, 0.75, by = 0.05),
      nrow(pars)
    ),
    actPropMiss = prop_miss
  )
  
}

# # find any problem cases and drop from both the parameter
# # dataframe and the data list
# probs <- unique(c(
#   which(
#     sapply(dat_flat, function(x){sum(is.infinite(x))}) > 0
#   ),
#   which(
#     sapply(dat_flat, function(x){sum(is.nan(x))}) > 0
#   ),
#   which(sapply(dat_flat, function(x){sum(x == 0, na.rm = T)}) > 0)
# ))

# dat_flat_complete <- dat_flat[-probs]
# pars_complete <- pars_full[-probs, ]

# time testing only
dat_flat <- dat_flat[as.numeric(in_args[5]):as.numeric(in_args[6])]
pars_full <- pars_full[as.numeric(in_args[5]):as.numeric(in_args[6]),]


### fitting the models ###

system.time({
  
# store results in a vector
results_vec <- numeric(8)
names(results_vec) <- c("estim_r","estim_alpha","se_r","se_alpha","lower_r","lower_alpha","upper_r","upper_alpha")

# only use MI

res1 <- lapply(
  X = dat_flat,
  FUN = fit_ricker_MI
)

res1



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

# add in forecasts here
forecast_rmse=function(ralpha,trueTS){
  r=ralpha[1]
  alpha=ralpha[2]
  l_ts=length(trueTS)
  # we assume that the 1st value in trueTS is the given starting value for all
  pred_TS=numeric(l_ts-1) 
  pred_TS[1]=trueTS[1]
  for(i in 1:(l_ts-1)){
    pred_TS[i+1]=pred_TS[i]*exp(r-alpha*pred_TS[i])
  }
  
  RMSE=sqrt(mean((pred_TS[2:l_ts] - trueTS[2:l_ts])^2))
  return(RMSE)
}

cur_index=as.numeric(in_args[5])

cur_sim=pars_full$SimNumber[1]
trueTS1=dat0[[cur_sim]]$y
trimmedlength=length(dat_flat[[1]])
trueTS=trueTS1[(trimmedlength+1):length(trueTS1)]

if(is.na(results_vec[1])|is.na(results_vec[2])){
  forecasts_MI=NA
} else {
  forecasts_MI=forecast_rmse(c(results_vec[1],results_vec[2]),trueTS=trueTS)
}

results_vec[9]=forecasts_MI
names(results_vec)[9]="forecasts_MI"



if(is.na(in_args[8])){
  ran20=F
} else {
  ran20=as.numeric(in_args[8])>18
}
if(is.na(results_vec[2])&results_vec[1]!="population extinction"&!ran20){
  # don't save if result was NA, and reason wasn't extinction and we haven't run 20 yet
} else {
  # save results to file
  write.csv(cbind(pars_full,as.data.frame(t(results_vec))),file = here(in_args[4]))
}


