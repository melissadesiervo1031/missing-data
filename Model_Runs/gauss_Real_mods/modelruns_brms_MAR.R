# Load packages ## 
#make sure these are already in the folder on supercomputer where I need them ##

#.libPaths("/pfs/tc1/home/astears/R/x86_64-pc-linux-gnu-library/4.2")

library(tidyverse)
library(brms)

# This script will run 3 ARIMA functions (drop missing, Kalman, Multiple imputations 
#over a nested list with increasing prop missing, over 1000+ simulations ###)

#CurSim = like a loop ##
# 
# CurSim <- commandArgs(trailingOnly = TRUE) #Look at command line arguments only after the R script
# CurSim <- as.numeric(CurSim)
# CurSim <- CurSim + 1 # since the Slurm array is 0 indexed

## read in the data

gauss_auSable_randMiss <- readRDS(here("data/missingDatasets/gauss_real_auSable_randMiss.rds"))
gauss_badger_randMiss <- readRDS(here("data/missingDatasets/gauss_real_badger_randMiss.rds"))
au_sable_river_full <- read_csv(here("data/au_sable_river_prepped.csv"))
badger_mill_creek_full <- read_csv(here("data/badger_mill_Creek_prepped.csv"))

# make file for output beforehand in supercomputer folder 
# will put them all together after all run, using the command line
#OutFile <- paste0("gauss_real/gauss_real_MAR_brms_modelResults_normPriorNB/", CurSim, "brmsvals.csv")
#OutFile_preds <- paste0("gauss_real/gauss_real_MAR_brms_modelResults_normPriorNB/", CurSim, "brmspreds.csv")

#########################################################################################
### MY ARIMA FUNCTIONS #####
##########################################################################################
### Function to fit a BRMS model on a time series ###
fit_brms_model <- function(sim_list, sim_pars, 
                           iter = 4000, include_missing = FALSE,
                           forecast = TRUE, forecast_days = 365,
                           dat_full ){
  simmissingdf <-lapply(X = sim_list, 
                        FUN = function(X) cbind.data.frame(GPP = X, 
                                                           light = sim_pars$light, 
                                                           discharge = sim_pars$Q))
  
  if(forecast){
    simmissingdf <- lapply(simmissingdf, function(df) {
      df[1:(nrow(df)-forecast_days), ]  # Remove to save these for forecasting
    })
  }
  
  # Make the model formula and priors
  bform <- brms::bf(GPP | mi() ~ light + Q + ar(p = 1))
  bprior <- c(prior(normal(0,1), class = 'ar'),
              prior(normal(0,5), class = 'b'))
  
  # fit model to list of datasets
  bmod <- brms::brm_multiple(bform, data = simmissingdf, 
                             prior = bprior, iter = iter, 
                             combine = FALSE)
  
  extract_brms_pars <- function(bfit, include_missing = FALSE){
    bsum <- brms::posterior_summary(bfit, probs = c(0.025, 0.5, 0.975))
    bsum <- as.data.frame(bsum) %>%
      mutate(parameter = row.names(bsum)) %>%
      filter(!(parameter %in% c('lprior', 'lp__'))) %>%
      mutate(parameter = case_when(parameter == 'ar[1]' ~ 'phi',
                                   TRUE ~ parameter)) %>%
      select(parameter, mean = Estimate, sd = Est.Error, 
             '2.5%' = Q2.5, '50%' = Q50, '97.5%' = Q97.5)
    
    if(!include_missing){
      bsum <- bsum[grep('^Ymi', bsum$parameter, invert = TRUE),]
    }
    
    row.names(bsum) <- NULL
    
    return(bsum)
  }
  
  bpars <- lapply(bmod, extract_brms_pars, include_missing = include_missing)
  names(bpars) <- names(simmissingdf)
  
  if(forecast){  
    dat_forecast <- dat_full %>%
      slice((nrow(dat_full)-forecast_days):nrow(dat_full)) %>%
      select(date, GPP, light, Q)
    
    predictions <- lapply(bmod, function(mod){
      predict(mod, newdata = dat_forecast[,-2]) %>%
        as.data.frame() %>% mutate(date = dat_forecast$date,
                                   GPP = dat_forecast$GPP)
    })
    
    names(predictions) <- names(simmissingdf)
    return(list(brms_forecast = predictions,
                brms_pars = bpars,
                sim_params = sim_pars))
  }
  
  return(list(brms_pars = bpars,
              sim_params = sim_pars))
  
}


# Fit models for au sable river data ---------------------------------------------
dirname <- c("./data/model_results/gauss_real_MAR_brms_modResults/badgerMill/")
for (i in seq_along(gauss_auSable_randMiss)) {
  CurSim <- i
  OutFile_params <- paste(dirname, CurSim, "brmsvals.csv", sep = "")
  OutFile_preds <- paste(dirname, CurSim, "brmspreds.csv", sep = "")
  
  # make sure sim_list and sim_params are not pointing to the whole list, which starts w character elements
  sim_list<- gauss_auSable_randMiss[[CurSim]]$y
  nms <- stringr::str_which(names(sim_list), "^prop")
  sim_list <- sim_list[nms]
  
  brms_MAR <- fit_brms_model(sim_list = sim_list,
                             sim_pars = gauss_auSable_randMiss[[CurSim]]$sim_params,
                             forecast = TRUE, forecast_days = 365, 
                             dat_full = au_sable_river_full)
  
  brms_MAR_df <- map_df(brms_MAR$brms_pars, ~as.data.frame(.x),
                        .id = "missingprop_autocor")
  brms_MAR_df$missingness <- 'MAR'
  brms_MAR_df$type <- 'brms'
  brms_MAR_df$run_no <- CurSim
  
  brms_MAR_preds <- map_df(brms_MAR$brms_forecast, ~as.data.frame(.x),
                           .id = "missingprop_autocor")
  brms_MAR_preds$missingness <- 'MAR'
  brms_MAR_preds$type <- 'brms'
  brms_MAR_preds$run_no <- CurSim
  #################################################
  # Write the output to the folder which will contain all output files as separate csv
  #    files with a single line of data.
  write_csv(brms_MAR_df, file = OutFile_params)
  write_csv(brms_MAR_preds, file = OutFile_preds)
}

# Fit models for badger mill creek data -----------------------------------
dirname <- c("./data/model_results/gauss_real_MAR_brms_modResults/badgerMill")
for (i in seq_along(gauss_badger_randMiss)) {
  CurSim <- i
  OutFile_params <- paste(dirname, CurSim, "brmsvals.csv", sep = "")
  OutFile_preds <- paste(dirname, CurSim, "brmspreds.csv", sep = "")
  
  # make sure sim_list and sim_params are not pointing to the whole list, which starts w character elements
  sim_list<- gauss_badger_randMiss[[CurSim]]$y
  nms <- stringr::str_which(names(sim_list), "^prop")
  sim_list <- sim_list[nms]
  
  brms_MAR <- fit_brms_model(sim_list = sim_list,
                             sim_pars = gauss_badger_randMiss[[CurSim]]$sim_params,
                             forecast = TRUE, forecast_days = 365, 
                             dat_full = badger_mill_creek_full)
  
  brms_MAR_df <- map_df(brms_MAR$brms_pars, ~as.data.frame(.x),
                        .id = "missingprop_autocor")
  brms_MAR_df$missingness <- 'MAR'
  brms_MAR_df$type <- 'brms'
  brms_MAR_df$run_no <- CurSim
  
  brms_MAR_preds <- map_df(brms_MAR$brms_forecast, ~as.data.frame(.x),
                           .id = "missingprop_autocor")
  brms_MAR_preds$missingness <- 'MAR'
  brms_MAR_preds$type <- 'brms'
  brms_MAR_preds$run_no <- CurSim
  #################################################
  # Write the output to the folder which will contain all output files as separate csv
  #    files with a single line of data.
  write_csv(brms_MAR_df, file = OutFile_params)
  write_csv(brms_MAR_preds, file = OutFile_preds)
}



# compile outputs ---------------------------------------------------------
# compile for au sable river
# compile script output and save 
dirname <- "./data/model_results/gauss_real_MAR_brms_modResults/auSable/"
predNames <- list.files(dirname, pattern = "preds.csv")
predsAll <- map_df(predNames, function(x) {
  read_csv(paste0(dirname, x))
})
#write.csv(predsAll, file = "./data/model_results/gauss_real_MAR_arima_FORECASTpreds.csv")
write.csv(predsAll, file = paste(dirname, "gauss_auSable_real_MAR_brms_FORECASTpreds.csv", sep = ""))

valNames <- list.files(dirname, pattern = "vals.csv")
valsAll <- map_df(valNames, function(x) {
  read_csv(paste0(dirname, x))
})
#write.csv(valsAll, file = "./data/model_results/gauss_real_MAR_arima_FORECASTvals.csv")
write.csv(valsAll, file = paste(dirname, "gauss_auSable_real_MAR_arima_FORECASTvals.csv", sep = ""))


# compile for badger mill creek
# compile script output and save 
dirname <- "./data/model_results/gauss_real_MAR_brms_modResults/badgerMill/"
predNames <- list.files(dirname, pattern = "preds.csv")
predsAll <- map_df(predNames, function(x) {
  read_csv(paste0(dirname, x))
})
#write.csv(predsAll, file = "./data/model_results/gauss_real_MAR_arima_FORECASTpreds.csv")
write.csv(predsAll, file = paste(dirname, "gauss_badger_real_MAR_brms_FORECASTpreds.csv", sep = ""))

valNames <- list.files(dirname, pattern = "vals.csv")
valsAll <- map_df(valNames, function(x) {
  read_csv(paste0(dirname, x))
})
#write.csv(valsAll, file = "./data/model_results/gauss_real_MAR_arima_FORECASTvals.csv")
write.csv(valsAll, file = paste(dirname, "gauss_badger_real_MAR_arima_FORECASTvals.csv", sep = ""))
