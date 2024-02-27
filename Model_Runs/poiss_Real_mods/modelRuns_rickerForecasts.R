#///////////////////////////////////////////////////
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
#///////////////////////////////////////////////////

library(here)
library(parallel)
library(tidyverse)

# source functions
f_list <- list.files("./Functions/", full.names = T)
#f_list <- list.files(here("Functions/"), full.names = T)
lapply(f_list, source)


# Prepping data -----------------------------------------------------------

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
  temp <- x[1:54]
  # remove NAs at the end of the dataset (not a realistic scenario... will this
  # be problematic? because the lenght of time series will be different?)
  if (sum(is.na(temp)) > 0) {
    if (max(which(!is.na(temp))) < 54) {
      temp2 <- temp[1:max(which(!is.na(temp)))]
    } else {
      temp2 <- temp
    }
    return(temp2)
  } else {
     return(temp) 
    }
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
  
# calculate the length of the trimmed time series
namesDF$n_trimmed <- unlist(map(dat_trimmed, function(x) length(x)))

# Fitting Models ----------------------------------------------------------
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

# get names of DA_fits data
DA_fits_df <- data.frame("name" = names(DA_fits),
                         "DA_fits" = NA)
DA_fits_df$DA_fits <- DA_fits
#drop any duplicated rows
DA_fits_df <- DA_fits_df[!(DA_fits_df$name %in% dups$oldName),]

# add all values together
allDat <- full_join(dropNA_fits_df, dropNAcc_fits_df, by = c("name")) %>% 
 left_join(MI_fits_df) %>% 
  left_join(EM_fits_df) %>% 
  left_join(DA_fits_df)
# get new names w/ accurate autocorr and propMiss data
allDat <- allDat %>% 
  rename(oldName = name) %>% 
  left_join(namesDF, by = "oldName") %>% 
  select(newName_trimmed, oldName, actAutoCorr_trim, actPropMiss_trim, drop_fits, cc_fits, MI_fits, 
         EM_fits, DA_fits) %>% 
  rename(newName = newName_trimmed)
# add in trimmed time series
names(dat_trimmed) <- NULL
allDat$trimmed_ts <- dat_trimmed

# Predict remaining values ------------------------------------------------
# use 263 as the starting value, even though some of the missing datasets end in an NA
#forecast_outputs <- 
allDat$forecasts <- (apply(allDat, MARGIN = 1, FUN = function(x) {
  # get starting value (last value of input time series)
  N_tminus1 <-  x$trimmed_ts[length(x$trimmed_ts)]
  counter <- length(x$trimmed_ts)
  # make empty matrix to store values
   outDat <- matrix(nrow = 59, ncol = 6)
   # define starting values for each value
  N_now_dropNA <- N_now_cc <- N_now_MI <- N_now_EM <-  N_now_DA <- N_tminus1
  outDat[counter,] <- c(counter, N_now_dropNA, N_now_cc, N_now_MI, N_now_EM, N_now_DA)
  while(counter < 59) {
    # define the "next value" for each model type
    N_next_dropNA <- rpois(1, lambda = N_now_dropNA * exp(x$drop_fits[[1]]["r"] - x$drop_fits[[1]]["alpha"]*N_now_dropNA))
    N_next_cc <- rpois(1, lambda = N_now_cc * exp(x$drop_fits[[1]]["r"] - x$drop_fits[[1]]["alpha"]*N_now_cc))
    N_next_MI <- rpois(1, lambda = N_now_MI * exp(x$drop_fits[[1]]["r"] - x$drop_fits[[1]]["alpha"]*N_now_MI))
    N_next_EM <- rpois(1, lambda = N_now_EM * exp(x$drop_fits[[1]]["r"] - x$drop_fits[[1]]["alpha"]*N_now_EM))
    N_next_DA <- rpois(1, lambda = N_now_DA * exp(x$drop_fits[[1]]["r"] - x$drop_fits[[1]]["alpha"]*N_now_DA))
    # save values w/ the time step
    outDat[counter+1,] <- c(counter+1, N_next_dropNA, N_next_cc, N_next_MI, N_next_EM, N_next_DA)
    # update "now" values
    N_now_dropNA <- N_next_dropNA
    N_now_cc <- N_next_cc
    N_now_MI <- N_next_MI
    N_now_EM <- N_next_EM
    N_now_DA <- N_next_DA
    #udpate counter
    counter <- counter+1
  }
  # add a column to "outDat" that has the intput time series
  
  outDat <- cbind(outDat, 
                  # add in trimmed ds
                  c(x$trimmed_ts,
                            rep(NA, length.out = (61- length(x$trimmed_ts[[1]]))),recursive = TRUE),
                  # add in complete time series
                  c(dat$pois_real_randMiss_autoCor_0$y$y, NA, recursive = TRUE)
                  ) 
  outDat <- as.data.frame(outDat) 
  names(outDat) <- c("timeStep", "dropNA_est", "dropCC_est", "MI_est", "EM_est", "DA_est", "trimmedInput_ts", "real_ts")
  
  return(outDat)
  counter <- NULL
  outDat <- NULL
}
 ))

#plot the results for one output 
plot(x = as.numeric(row.names(allDat[1,"forecasts"][[1]])), 
     y = allDat[1,"forecasts"][[1]][,"real_ts"], 
     type = "l")
lines(x = 1:59, 
      y = allDat[1,"forecasts"][[1]][,"dropNA_est"], col = "red")
lines(x = 1:59, 
      y = allDat[1,"forecasts"][[1]][,"dropCC_est"], col = "orange")
lines(x = 1:59, 
      y = allDat[1,"forecasts"][[1]][,"MI_est"], col = "green")
lines(x = 1:59, 
      y = allDat[1,"forecasts"][[1]][,"EM_est"], col = "blue")
lines(x = 1:59, 
      y = allDat[1,"forecasts"][[1]][,"DA_est"], col = "purple")


# Calculate RMSE  ---------------------------------------------------------
allDat$RMSE <- lapply(allDat$forecasts, FUN = function(x) {
  # get starting index
  startTimeStep <- which(!is.na(x$timeStep))[2]
  # calculate RMSE for dropNA
  RMSE_dropNA <- sqrt(mean((x$dropNA_est[startTimeStep:59] - x$real_ts[startTimeStep:59])^2))
  # calculate RMSE for dropNA_cc
  RMSE_dropCC <- sqrt(mean((x$dropCC_est[startTimeStep:59] - x$real_ts[startTimeStep:59])^2))
  # calculate RMSE for MI
  RMSE_MI <- sqrt(mean((x$MI_est[startTimeStep:59] - x$real_ts[startTimeStep:59])^2))
  # calculate RMSE for EM
  RMSE_EM <- sqrt(mean((x$EM_est[startTimeStep:59] - x$real_ts[startTimeStep:59])^2))
  # calculate RMSE for DA
  RMSE_DA <- sqrt(mean((x$DA_est[startTimeStep:59] - x$real_ts[startTimeStep:59])^2))
  return(c("dropNA" = RMSE_dropNA, "dropCC" = RMSE_dropCC, "MI" = RMSE_MI, "EM" = RMSE_EM, "DA" = RMSE_DA))
  }
)


## save the output     
saveRDS(allDat, file = "./data/model_results/RickerForecast_resultTableAll.rds")
