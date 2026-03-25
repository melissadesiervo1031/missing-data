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

# REVISION NEEDS: 
# 1) Add in forecasting
# NOTES: how did we do forecasting for the real data? Maybe run separate script- made into useful for this?
# 2) why don't the numbers match up for the MinMaxMiss data set- should be less because no autocorrelation?
# ANS: because we also include the full data set as an entry, meaning 16 instead of 15 entries per simulation
# 3) Compare to model run on full data set with 0% missing
# ANS: we include the full in MinMaxMiss, so we should be good if we use that to forecast but not to do the missingness- see if 0 missing can be skipped for missingness trials
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
#in_args=c("data/missingDatasets/pois_sim_MinMaxMiss.rds", "data/missingDatasets/pois_sim_params.rds", 2, "Model_Runs/Ricker_example1.rds", 15001, 15007, "drop")
#in_args=c("data/missingDatasets/pois_sim_randMiss_B.rds", "data/missingDatasets/pois_sim_params.rds", 2, "Model_Runs/Ricker_example1.rds", 1, 5)

in_args <- commandArgs(trailingOnly = T)
cat(in_args)
# read in datafile
dat <- readRDS(here(in_args[1]))
pars <- readRDS(here(in_args[2]))
dat0=readRDS(here("data/ricker_0miss_datasets.rds"))
cat("data has been loaded")

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
      seq(0.00, 0.75, by = 0.05),# RandMiss includes the full data set, so include 0 here
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
      seq(0.00, 0.75, by = 0.05), # MinMaxMiss includes the full data set, so include 0 here
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
dat_flat_o=dat_flat
dat_flat <- dat_flat[as.numeric(in_args[5]):as.numeric(in_args[6])]
pars_full <- pars_full[as.numeric(in_args[5]):as.numeric(in_args[6]),]


### fitting the models ###

# define methods to be used
if(is.na(in_args[7])){
  methods <- paste0(
    c("drop", "cc", "EM","DA")
  )
} 



system.time({

# make cluster for parallel computing
cl <- parallelly::makeClusterPSOCK(as.numeric(in_args[3]))

clusterEvalQ(cl = cl, expr = {
  library(here)
  f_list <- list.files(here("Functions/"), full.names = T)
  f_list=f_list[-grep("README",f_list)]
  lapply(f_list, source)
})

# store results in a list
results_list <- vector(mode = "list", length = length(methods))
names(results_list) <- paste0(
  methods, "_fits"
)



for(i in 1:length(methods)){
  fit_fun <- eval(parse(text = paste0("fit_ricker_", methods[i])))
  
  clusterExport(cl, "fit_fun", envir = environment())
  
  safe_fun <- function(x) {
    if (all(is.na(x))) {
      return(NA)
    }
    fit_fun(x)
  }
  results_list[[i]] <- parLapply(
    cl = cl,
    X = dat_flat,
    fun = safe_fun
  )
}


stopCluster(cl)


# compile results into a tibble
results <- cbind(
  as_tibble(pars_full),
  as_tibble(results_list)
)
})

# forecasts

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

lengths=numeric(10)
for(i in 1:10){
  lengths[i]=length(dat_flat_o[[i]])
}
u_lengths=unique(lengths)
frequencies=numeric(length(u_lengths))
for(i in 1:length(u_lengths)){
  frequencies[i]=length(which(lengths==u_lengths[i]))
}
trimmedlength=unique(lengths)[which.max(frequencies)]

forecasts=matrix(data=NA,nrow=nrow(results),ncol=length(methods))

for(i in 1:nrow(results)){

  for(j in 1:length(methods)){
    column_cur=grep(paste0(methods[j],"_fits"),colnames(results))
    ralpha=results[i,column_cur][[1]][[1]] 
    cur_sim=results[i,1]
    trueTS1=dat0[[cur_sim]]$y
    trueTS=trueTS1[(trimmedlength+1):length(trueTS1)]
    forecasts[i,j]=forecast_rmse(ralpha=ralpha,trueTS=trueTS)
  }
  
}


results <- cbind(
  results,
  forecasts
)

colnames(results)[(ncol(results)-length(methods)+1):ncol(results)]=paste0("forecast_RMSE_",methods)


# save results to file
saveRDS(results, file = here(in_args[4]))


