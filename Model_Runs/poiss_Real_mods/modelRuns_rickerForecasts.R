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
dat <- readRDS("./data/missingDatasets/pois_real_randMiss.rds")


# count number of missingness proportions (within each level of autocorrelation)
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


### leaving out data for forecasting (last 5 years)
dat_complete <- dat_flat
dat_trimmed <- map(dat_flat, function(x) {
  x[1:54]
})

# double check the autocorrelation of each timeseries
autoCorr_actual_full <- map(dat_flat, 
                       function(x) {
                         x[!is.na(x)] <- 0
                         x[is.na(x)] <- 1 
                         return(round(acf(x, plot = FALSE)$acf[2],2))
                       }
)
propMissing_actual_full <- map(dat_flat,
                          function(x){
                            return(round(sum(is.na(x))/length(x),2))
                          })

# make a lookup table w/ the new names (to fix the names later)
namesDF <- data.frame("oldName" = names(dat_flat),
                      "actAutoCorr_full" = unlist(autoCorr_actual_full),
                      "actPropMiss_full" = unlist(propMissing_actual_full))
# make new names
namesDF$newName_partial <- apply(as.matrix(str_split(namesDF$oldName, pattern = "_", simplify = TRUE)[,1:3]), MARGIN = 1, 
                                 function(x) {
                                   str_flatten(x, collapse = "_")
                                 })
namesDF$newName_full <- paste0(namesDF$newName_partial, 
                          "autoCorr_", namesDF$actAutoCorr_full, 
                          "_propMiss_", namesDF$actPropMiss_full)

# double check the autocorrelation of each timeseries
autoCorr_actual_trim <- map(dat_trimmed, 
                            function(x) {
                              x[!is.na(x)] <- 0
                              x[is.na(x)] <- 1 
                              return(round(acf(x, plot = FALSE)$acf[2],2))
                            }
)
propMissing_actual_trim <- map(dat_trimmed,
                               function(x){
                                 return(round(sum(is.na(x))/length(x),2))
                               })

# make a lookup table w/ the new names (to fix the names later)
namesDF <- left_join(namesDF, 
                          data.frame("oldName" = names(dat_trimmed),
                                     "actAutoCorr_trim" = unlist(autoCorr_actual_trim),
                                     "actPropMiss_trim" = unlist(propMissing_actual_trim)), 
                          by = "oldName")
namesDF$newName_trimmed <- paste0(namesDF$newName_partial, 
                                    "autoCorr_", namesDF$actAutoCorr_trim, 
                                    "_propMiss_", namesDF$actPropMiss_trim)
  
### fitting the models ###
# few enough that I think I can run locally... 

# dropNA simple case
dropNA_fits <- map(dat_trimmed, function(x) fit_ricker_drop(x, fam = "poisson"))

# dropNA complete case
dropNAcc_fits <- map(dat_flat, function(x) fit_ricker_cc(y = x, fam = "poisson")) 

# multiple imputations
MI_fits <- map(dat_flat, function(x) fit_ricker_MI(y = x))

# expectation maximization
EM_fits <- map(dat_flat, function(x) fit_ricker_EM(y = x))

# data augmentation # takes too long to run locally... 
DA_fits <- map(dat_flat, function(x) fit_ricker_DA(y = x))

## convert to a data.frame, and add together
# get the names of any duplicated input names
dups <- namesDF[namesDF$oldName %in% namesDF[duplicated(namesDF$oldName),"oldName"],]

# get names of dropNA_fits data
dropNA_fits_df <- data.frame("name" = names(dropNA_fits),
                             "drop_fits" = NA)
dropNA_fits_df$drop_fits <- dropNA_fits
#drop any duplicated rows
dropNA_fits_df <- dropNA_fits_df[!(dropNA_fits_df$name %in% dups$oldName),]

# get names of dropNAcc_fits data
dropNAcc_fits_df <- data.frame("name" = names(dropNAcc_fits),
                               "cc_fits" = NA)
dropNAcc_fits_df$cc_fits <- dropNAcc_fits
#drop any duplicated rows
dropNAcc_fits_df <- dropNAcc_fits_df[!(dropNAcc_fits_df$name %in% dups$oldName),]

# get names of MI_fits data
MI_fits_df <- data.frame("name" = names(MI_fits),
                         "MI_fits" = NA)
MI_fits_df$MI_fits <- MI_fits
#drop any duplicated rows
MI_fits_df <- MI_fits_df[!(MI_fits_df$name %in% dups$oldName),]

# get names of EM_fits data
EM_fits_df <- data.frame("name" = names(EM_fits),
                         "EM_fits" = NA)
EM_fits_df$EM_fits <- EM_fits
#drop any duplicated rows
EM_fits_df <- EM_fits_df[!(EM_fits_df$name %in% dups$oldName),]

# # get names of DA_fits data
# DA_fits_df <- data.frame("name" = names(DA_fits),
#                          "DA_fits" = NA)
# DA_fits_df$DA_fits <- DA_fits
# #drop any duplicated rows
# DA_fits_df <- DA_fits_df[!(DA_fits_df$name %in% dups$oldName),]

# add all values together
allDat <- full_join(dropNA_fits_df, dropNAcc_fits_df, by = c("name")) %>% 
  left_join(MI_fits_df) %>% 
  left_join(EM_fits_df)
# get new names w/ accurate autocorr and propMiss data
allDat <- allDat %>% 
  rename(oldName = name) %>% 
  left_join(namesDF, by = "oldName") %>% 
  select(newName, actAutoCorr, actPropMiss, drop_fits, cc_fits, MI_fits, EM_fits) %>% 
  rename(name = newName)
# get simulation number in it's own column
allDat <- allDat %>% 
  mutate(simNumber = as.integer(str_split(allDat$name, "_", simplify = TRUE)[,2] %>% 
                                  str_extract("\\d+"))) %>% 
  select(name, simNumber, actAutoCorr, actPropMiss, drop_fits, cc_fits, MI_fits, EM_fits)

# add back in the simulation parameters 
allDat <- allDat %>% 
  left_join(pars, by = c("simNumber" = "SimNumber")) 


## save the output     
saveRDS(allDat, file = "./data/model_results/RickerExtinct_resultTableAll.rds")
