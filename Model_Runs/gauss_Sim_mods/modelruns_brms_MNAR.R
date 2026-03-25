# Load packages ## 
#make sure these are already in the folder on supercomputer where I need them ##

#.libPaths("/pfs/tc1/home/astears/R/x86_64-pc-linux-gnu-library/4.2")

library(tidyverse)
library(brms)

# This script will run data augmentation models using the BRMS package 
#over a nested list with increasing prop missing, over 1000+ simulations ###)

#CurSim = like a loop ##

# CurSim <- commandArgs(trailingOnly = TRUE) #Look at command line arguments only after the R script
# CurSim <- as.numeric(CurSim)
# CurSim <- CurSim + 1 # since the Slurm array is 0 indexed

## read in the autocor_01 list ##

gauss_sim_MNAR <- readRDS("data/missingDatasets/gauss_sim_minMaxMiss.rds")

# make file for output beforehand in supercomputer folder 
# will put them all together after all run, using the command line
OutFile <- paste0("./data/model_results/gauss_sim_minMax_modelResults/")

# OR create directory for local output  - COMMENT OUT FOR HPC!
# dir.create("./data/model_results/gauss_sim_minMax_modelResults/", recursive = TRUE, showWarnings = FALSE)


#########################################################################################
### MY BRMS FUNCTIONS #####
##########################################################################################

fit_brms_model <- function(sim_list, sim_pars, 
                           iter = 4000, include_missing = FALSE,
                           forecast = TRUE, forecast_days = 73,
                           dat_full){
  
  simmissingdf <-lapply(X = sim_list, 
                        FUN = function(X) cbind.data.frame(GPP = X[1:292], 
                                                           light = sim_pars$X[,2][1:292], 
                                                           discharge = sim_pars$X[,3][1:292]))
  
  
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
      slice((nrow(dat_full)-forecast_days):nrow(dat_full)) %>%  # Minor point (CT): this gives 74 forecast days (obs 292-365 inclusive); for 73 forecast days, would need to change to (nrow(dat_full) - forecast_days + 1)
      select(date, GPP, light, discharge)
    
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
#### MODEL RUN BRMS DROP ##############
#########################################################
# slightly increase the control of size of object fitted by parallel brms 
options(future.globals.maxSize = 1.0 * 1e10)

# fit models
for (i in 976:length(gauss_sim_MNAR)) {
  CurSim <- i

  ###
  # Commented out as they were duplicated below, under 'SAVE' - CT
  # OutFile_params <- paste(OutFile, CurSim, "brmsvals.csv", sep = "")
  # OutFile_preds <- paste(OutFile, CurSim, "brmspreds.csv", sep = "")
  
  # make sure sim_list and sim_params are not pointing to the whole list, which starts w character elements
  sim_list<- gauss_sim_MNAR[[CurSim]]$y
  nms <- stringr::str_which(names(sim_list), "^prop")
  sim_list <- sim_list[nms]
  
  # the 'full' dataset is the first list in the 'sim_list" 
  dat_full <- data.frame("date" = 1:365,
                         "GPP" = gauss_sim_MNAR[[CurSim]]$y$y_noMiss, 
                         "light" = gauss_sim_MNAR[[CurSim]]$sim_params$X[,2],  #changed from 1 to CurSim, assuming we don't want to use sim1 data for every iteration  - CT
                         "discharge" = gauss_sim_MNAR[[CurSim]]$sim_params$X[,3] #changed from 1 to CurSim, for above reasons - CT
  )

  # fit models
  brms_MNAR <- fit_brms_model(sim_list = sim_list,
                             sim_pars = gauss_sim_MNAR[[CurSim]]$sim_params,
                             forecast = TRUE, forecast_days = 73, 
                             dat_full = dat_full)

  ########### formatting for figure #############
  brms_MNAR_df <- map_df(brms_MNAR$brms_pars, ~as.data.frame(.x),
                        .id = "missingprop_autocor")
  brms_MNAR_df$missingness <- 'MNAR'
  brms_MNAR_df$type <- 'brms'
  brms_MNAR_df$run_no <- CurSim
  
  brms_MNAR_preds <- map_df(brms_MNAR$brms_forecast, ~as.data.frame(.x),
                           .id = "missingprop_autocor")
  brms_MNAR_preds$missingness <- 'MNAR'
  brms_MNAR_preds$type <- 'brms'
  brms_MNAR_preds$run_no <- CurSim
  ###################################################
  #### SAVE #########
  #################################################
  # Write the output to the folder which will contain all output files as separate csv
  #    files with a single line of data.
  OutFile_params <- paste(OutFile, CurSim, "brmsvals.csv", sep = "")
  OutFile_preds <- paste(OutFile, CurSim, "brmspreds.csv", sep = "")
 
  #################################################
  # Write the output to the folder which will contain all output files as separate csv
  #    files with a single line of data.
  write_csv(brms_MNAR_df, file = OutFile_params)
  write_csv(brms_MNAR_preds, file = OutFile_preds)
  
}

# fit the model w/ no missingness
for (i in 491:length(gauss_sim_MNAR)) {
  CurSim <- i
  
  ###
  # Commented out as they were duplicated below, under 'SAVE' - CT
  # OutFile_params <- paste(OutFile, CurSim, "brmsvals.csv", sep = "")
  # OutFile_preds <- paste(OutFile, CurSim, "brmspreds.csv", sep = "")
  
  # make sure sim_list and sim_params are not pointing to the whole list, which starts w character elements
  sim_list<- gauss_sim_MNAR[[CurSim]]$y
 # nms <- stringr::str_which(names(sim_list), "^prop")
  sim_list <- sim_list["y_noMiss"]
  
  # the 'full' dataset is the first list in the 'sim_list" 
  dat_full <- data.frame("date" = 1:365,
                         "GPP" = gauss_sim_MNAR[[CurSim]]$y$y_noMiss, 
                         "light" = gauss_sim_MNAR[[CurSim]]$sim_params$X[,2],  #changed from 1 to CurSim, assuming we don't want to use sim1 data for every iteration  - CT
                         "discharge" = gauss_sim_MNAR[[CurSim]]$sim_params$X[,3]    #changed from 1 to CurSim, for above reasons - CT
  )
  )
  
  # fit models
  brms_MNAR <- fit_brms_model(sim_list = sim_list,
                              sim_pars = gauss_sim_MNAR[[CurSim]]$sim_params,
                              forecast = TRUE, forecast_days = 73, 
                              dat_full = dat_full)
  
  ########### formatting for figure #############
  brms_MNAR_df <- map_df(brms_MNAR$brms_pars, ~as.data.frame(.x),
                         .id = "missingprop_autocor")
  brms_MNAR_df$missingness <- 'MNAR'
  brms_MNAR_df$type <- 'brms'
  brms_MNAR_df$run_no <- CurSim
  
  brms_MNAR_preds <- map_df(brms_MNAR$brms_forecast, ~as.data.frame(.x),
                            .id = "missingprop_autocor")
  brms_MNAR_preds$missingness <- 'MNAR'
  brms_MNAR_preds$type <- 'brms'
  brms_MNAR_preds$run_no <- CurSim
  ###################################################
  #### SAVE #########
  #################################################
  # Write the output to the folder which will contain all output files as separate csv
  #    files with a single line of data.
  OutFile_params <- paste(OutFile, CurSim, "brmsvals_noMiss.csv", sep = "")
  OutFile_preds <- paste(OutFile, CurSim, "brmspreds_noMiss.csv", sep = "")
  
  #################################################
  # Write the output to the folder which will contain all output files as separate csv
  #    files with a single line of data.
  write_csv(brms_MNAR_df, file = OutFile_params)
  write_csv(brms_MNAR_preds, file = OutFile_preds)
  
}

# Once the job finishes, you can use the following command from within the folder
#    containing all single line csv files to compile them into a single csv file:
#     awk '(NR == 1) || (FNR > 1)' *brmsvals.csv > AllParams_brms.csv
#     awk '(NR == 1) || (FNR > 1)' *brmspreds.csv > AllPreds_brms.csv
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



