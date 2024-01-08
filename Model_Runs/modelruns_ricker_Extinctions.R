###################################################################
# This script is designed to run on a cluster and takes
# arguments from the command line including:
#   1. The file storing the list of datasets.ß
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
library(tidyverse)

# source functions
f_list <- list.files(here("Functions/"), full.names = T)
lapply(f_list, source)

# arguments from the shell
# Should include (1) data file (2) parameter file (3) number of nodes for cluster (4) save file location
# (5) beginnning index (6) ending index (7) optional- model list
in_args <- commandArgs(trailingOnly = T)
cat(in_args)

# read in datafile
dat <- readRDS("./data/missingDatasets/pois_sim_randMiss_extinctions.rds")
pars <- readRDS("./data/missingDatasets/pois_sim_params.rds")


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


## get rid of time series that are fewer than 5 observations (arbitrary thing here?)
dat_flat <- keep(dat_flat, function(x) length(x) >= 5)
# if(!is.null(names(dat))){
# 
#   # get input autocorrelation from names
#   autocorrs <- str_extract(
#     names(dat),
#     "autoCorr_\\d+"
#   ) %>% str_extract(
#     .,
#     pattern = "\\d+"
#   ) %>% as.numeric() / 100
# 
#   # get actual autocorrelation from names
#   autocorr_act <- map(
#     dat,
#     ~ str_extract(
#       names(.x$y),
#       "Corr_0.\\d+"
#     ) %>% str_extract(
#       pattern = "0.\\d+"
#     ) %>% as.numeric()
#   ) %>% unlist()
# 
#   # get number of different autocorrelations
#   n_autocorrs <- length(unique(autocorrs))
# 
#   pars_full <- pars[rep(1:nrow(pars), n_autocorrs), ]
#   pars_full <- pars_full[rep(1:nrow(pars_full), each = nmiss_props), ]
# 
#   pars_full <- pars_full %>% mutate(
#     id = 1:length(dat_flat),
#     autoCorr = rep(autocorrs, each = 15),
#     propMiss = rep(
#       seq(0.05, 0.75, by = 0.05),
#       nrow(pars) * n_autocorrs
#     ),
#     actAutoCorr = autocorr_act,
#     actPropMiss = prop_miss
#   )
# 
# } else {
# 
#   pars_full <- pars[rep(1:nrow(pars), each = nmiss_props), ]
#   pars_full <- pars_full %>% mutate(
#     id = 1:length(dat_flat),
#     propMiss = rep(
#       seq(0.05, 0.75, by = 0.05),
#       nrow(pars)
#     ),
#     actPropMiss = prop_miss
#   )
# 
# }
# 
# # # find any problem cases and drop from both the parameter
# # # dataframe and the data list
# # probs <- unique(c(
# #   which(
# #     sapply(dat_flat, function(x){sum(is.infinite(x))}) > 0
# #   ),
# #   which(
# #     sapply(dat_flat, function(x){sum(is.nan(x))}) > 0
# #   ),
# #   which(sapply(dat_flat, function(x){sum(x == 0, na.rm = T)}) > 0)
# # ))
# 
# # dat_flat_complete <- dat_flat[-probs]
# # pars_complete <- pars_full[-probs, ]
# 
# # time testing only
# dat_flat <- dat_flat[as.numeric(in_args[5]):as.numeric(in_args[6])]
# pars_full <- pars_full[as.numeric(in_args[5]):as.numeric(in_args[6]),]
# 

### fitting the models ###
 # few enough that I think I can run locally... 

# dropNA simple case
dropNA_fits <- map(dat_flat, function(x) fit_ricker_drop(x, fam = "poisson"))

# dropNA complete case
dropNAcc_fits <- map(dat_flat, function(x) fit_ricker_cc(y = x, fam = "poisson")) 

# multiple imputations

# expectation maximization

# data augmentation