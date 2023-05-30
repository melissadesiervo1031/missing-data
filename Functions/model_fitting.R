################################################################################
# Functions to fit models to simulated datasets with varying degrees of missing data
#

library(tidyverse)
library(brms)

################################################################################
# Data Augmentation in Stan using BRMS ####

#' fit_brms_model: use brms to fit a missing data model to simulated data
#'
#' This function fits a brms missing-data AR1 model to a simulated data set 
#' multiple times, once for each different level of missingness in the data.
#' It returns a list of parameters from each model fit as well as the parameters
#' used to simulate the data.
#'
#' @param sim_list a list of datasets simulated from the same set of parameters
#' each with a different level of missingness.
#' 
#' @param sim_pars a list of parameters that was used to simulate the data 
#' including the value for the AR1 coefficient (phi), the value of the covariate
#' effects (beta) and the matrix of simulated covariates (X)
#'
#' @param include_missing a boolean indicating if the parameter estimates for 
#' the missing data should be returned or not. Defaults to FALSE to match the 
#' arima function outputs.
#' 
#' @return A list of parameter estimates from fitting a stan AR1 model on each
#' simulated dataset.  
#'

fit_brms_model <- function(sim_list, sim_pars, include_missing = FALSE){
    simmissingdf <-lapply(X = sim_list, 
                          FUN = function(X) cbind.data.frame(GPP = X, 
                                                             light = sim_pars$X[,2], 
                                                             discharge = sim_pars$X[,3]))
    
    
    
    bform <- brms::bf(GPP | mi() ~ light + discharge + ar(p = 1))
    bmod <- brms::brm_multiple(bform, data = simmissingdf, combine = FALSE)
    
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
    
    return(list(brms_pars = bpars,
                sim_params = sim_pars))
    
}


# example code using this function:
# gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss.rds")
# GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]
# 
# brms_fit <- fit_brms_model(GPP_sim_MAR$y,
#                            GPP_sim_MAR$sim_params,
#                            include_missing = FALSE) 
