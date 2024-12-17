# Load packages ## 
# this file can run LOCALLY
library(tidyverse)
library(brms)

# This script will run data augmentation models using the BRMS package 
#over a nested list with increasing prop missing, over 1000+ simulations ###)

#CurSim = like a loop ##
CurSim <- 1 # since the Slurm array is 0 indexed

## read in the data ##
gauss_auSable_MNAR <- readRDS(here("data/missingDatasets/gauss_real_auSable_MinMaxMiss.rds"))
gauss_badger_MNAR <- readRDS(here("data/missingDatasets/gauss_real_badger_MinMaxMiss.rds"))
au_sable_river_full <- read_csv(here("data/au_sable_river_prepped.csv"))
badger_mill_creek_full <- read_csv(here("data/badger_mill_Creek_prepped.csv"))

# make file for output beforehand in supercomputer folder 
# will put them all together after all run, using the command line

#########################################################################################
### brms function #####
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
  bform <- brms::bf(GPP | mi() ~ light + discharge + ar(p = 1))
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
      select(date, GPP, light, Q) %>% 
      rename(discharge = Q)
    
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


#####################################################
#### MODEL run for au sable river##############
#########################################################
# Fit models for au sable river data ---------------------------------------------
dirname <- c("./data/model_results/gauss_real_MNAR_brms_modResults/auSable/")

OutFile_params <- paste(dirname,"brmsvals.csv", sep = "")
OutFile_preds <- paste(dirname, "brmspreds.csv", sep = "")

# make sure sim_list and sim_params are not pointing to the whole list, which starts w character elements
sim_list <- gauss_auSable_MNAR[[1]]$y
nms <- stringr::str_which(names(sim_list), "^prop")
sim_list <- sim_list[nms]

brms_MNAR <- fit_brms_model(sim_list = sim_list,
                            sim_pars = gauss_auSable_MNAR[[1]]$sim_params,
                            forecast = TRUE, forecast_days = 365, 
                            dat_full = au_sable_river_full)

brms_MNAR_df <- map_df(brms_MNAR$brms_pars, ~as.data.frame(.x),
                       .id = "missingprop_autocor")
brms_MNAR_df$missingness <- 'MNAR'
brms_MNAR_df$type <- 'brms'
brms_MNAR_df$run_no <- 1

brms_MNAR_preds <- map_df(brms_MNAR$brms_forecast, ~as.data.frame(.x),
                          .id = "missingprop_autocor")
brms_MNAR_preds$missingness <- 'MNAR'
brms_MNAR_preds$type <- 'brms'
brms_MNAR_preds$run_no <- 1
#################################################
# Write the output to the folder which will contain all output files as separate csv
#    files with a single line of data.
write_csv(brms_MNAR_df, file = OutFile_params)
write_csv(brms_MNAR_preds, file = OutFile_preds)


# Fit models for badger mill creek data -----------------------------------
dirname <- c("./data/model_results/gauss_real_MNAR_brms_modResults/badgerMill//")

OutFile_params <- paste(dirname,"brmsvals.csv", sep = "")
OutFile_preds <- paste(dirname, "brmspreds.csv", sep = "")

# make sure sim_list and sim_params are not pointing to the whole list, which starts w character elements
sim_list <- gauss_badger_MNAR[[1]]$y
nms <- stringr::str_which(names(sim_list), "^prop")
sim_list <- sim_list[nms]

brms_MNAR <- fit_brms_model(sim_list = sim_list,
                            sim_pars = gauss_badger_MNAR[[1]]$sim_params,
                            forecast = TRUE, forecast_days = 365, 
                            dat_full = badger_mill_creek_full)

brms_MNAR_df <- map_df(brms_MNAR$brms_pars, ~as.data.frame(.x),
                       .id = "missingprop_autocor")
brms_MNAR_df$missingness <- 'MNAR'
brms_MNAR_df$type <- 'brms'
brms_MNAR_df$run_no <- 1

brms_MNAR_preds <- map_df(brms_MNAR$brms_forecast, ~as.data.frame(.x),
                          .id = "missingprop_autocor")
brms_MNAR_preds$missingness <- 'MNAR'
brms_MNAR_preds$type <- 'brms'
brms_MNAR_preds$run_no <- 1
#################################################
# Write the output to the folder which will contain all output files as separate csv
#    files with a single line of data.
write_csv(brms_MNAR_df, file = OutFile_params)
write_csv(brms_MNAR_preds, file = OutFile_preds)



