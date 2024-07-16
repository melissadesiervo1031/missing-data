# Load packages ## 
# this file can run LOCALLY
library(tidyverse)
library(brms)

#over a nested list with increasing prop missing, over 1000+ simulations ###)

#CurSim = like a loop ##
CurSim <- 1 # since the Slurm array is 0 indexed

## read in the data ##
gauss_real_MNAR<- readRDS("./data/missingDatasets/gauss_real_MinMaxMiss.rds")
au_sable_full <- read_csv('./data/au_sable_river_prepped.csv')

# make file for output beforehand in supercomputer folder 
# will put them all together after all run, using the command line
OutFile <- ("gauss_real_MNAR_brms_FORECASTresults_normPriorNB.csv")
OutFile_preds <- ("gauss_real_MNAR_brms_FORECASTpreds_normPriorNB.csv")

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
      predict(mod, newdata = dat_forecast[,-2], n.ahead = forecast_days+1) %>%
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
#### MODEL RUN ARIMA DROP ##############
#########################################################

brms_MNAR <- fit_brms_model(sim_list = gauss_real_MNAR[[CurSim]]$y,
                           sim_pars = gauss_real_MNAR[[CurSim]]$sim_params,
                           forecast = TRUE, forecast_days = 365,
                           dat_full = au_sable_full)


#####

########### formatting for figure #############

brms_MNAR_df <- map_df(brms_MNAR$brms_pars, ~as.data.frame(.x),
                      .id = "missingprop_autocor")
brms_MNAR_df$missingness <- 'MAR'
brms_MNAR_df$type <- 'brms'
brms_MNAR_df$run_no <- CurSim

brms_MNAR_preds <- map_df(brms_MNAR$brms_forecast, ~as.data.frame(.x),
                         .id = "missingprop_autocor")
brms_MNAR_preds$missingness <- 'MAR'
brms_MNAR_preds$type <- 'brms'
brms_MNAR_preds$run_no <- CurSim

# fix "missingprop_autocor" column in brms_MNAR_preds so it shows the actual missingness quantity
unique(brms_MNAR_preds$missingprop_autocor)

namesActual <- data.frame("missingprop_autocor" = as.character(c(1:16)),
                          "actualName" =unique(names(gauss_real_MNAR[[1]]$y)))

brms_MNAR_preds <- brms_MNAR_preds %>% 
  #left_join(brms_MNAR_preds, namesActual) %>% 
  #select(-missingprop_autocor) %>% 
  #rename(missingprop_autocor = actualName) %>% 
  select("missingprop_autocor", "Estimate", "Est.Error", "Q2.5", "Q97.5", "date", "GPP", "missingness", "type", "run_no")
  
###################################################
#### SAVE #########
#################################################
# Write the output to the folder which will contain all output files as separate csv
#    files with a single line of data.
write_csv(brms_MNAR_df, file = paste0("./data/model_results/", OutFile))
write_csv(brms_MNAR_preds, file = paste0("./data/model_results/", OutFile_preds))


# Once the job finishes, you can use the following command from within the folder
#    containing all single line csv files to compile them into a single csv file:
#     awk '(NR == 1) || (FNR > 1)' *vals.csv > AllResults.csv
# The * is a wildcard character so the input to this will match any file within
#    your current folder that ends with vals.csv regardless of the rest of the filename.
#    These will then all be combined into a single file (AllResults.csv). The order
#    will be based on how linux orders the file names within the directory, so it 
#    might not match the original order of your parameter input file, but all the
#    entries will be there and it can be sorted later. Alternatively, you can name
#    your output files in a way in which the proper order will be enforced (e.g.,
#    if you will have a total of 100 jobs, you can name them all with 3 digits like
#    001_vals.csv, 002_vals.csv, etc.)
# Once you have combined all the single line csv files into your master results file,
#    you can remove them using the wildcard character again (e.g., rm *vals.csv)



