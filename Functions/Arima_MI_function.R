# Load packages
library(here)
library(stats)
library(forecast)
library(tidyverse)
library(lubridate)
library(Amelia)


#' Function that will fit the model using ARIMA while dealing with missing 
#' values using imputation via the Amelia package.
#'
#' @param sim_list a list of GPP data from multiple simulations
#' @param sim_pars a list of simulation parameters corresponding to the GPP simulated values
#' @param imputationsnum  the number of imputed datasets to create. This corresponds to the
#' \code{"m"} argument to the \code{"amelia"} function from the \code{"Amelia"} package.
#'
#' @return List of the ARIMA parameter estimates, errors, and true values from the simulations.
#'
#' @examples
#' 
#' gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_0.rds")
#' GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]
#' 
#' arima_MI <- fit_arima_MI(GPP_sim_MAR$y,GPP_sim_MAR$sim_params, imputationsnum=5)
#' 
fit_arima_MI <- function(sim_list, sim_pars, imputationsnum){
  
  # Extract the indices for the time series (1 through the maximum value) to pass to the 
  #   amelia function for imputation
  days <- seq(1, length(sim_list[[1]]))
  
  # For each simulation, combine the GPP values with the day indices and covariates in a single data frame
  sim_missing_df <- lapply(X = sim_list, 
                        FUN = function(X) cbind.data.frame(days= days,
                                                           GPP = X, 
                                                           light = sim_pars$X[,2], 
                                                           discharge = sim_pars$X[,3]))
  
  
  full_amelia <- lapply(X = sim_missing_df  , FUN = function(X)   amelia(X, ts="days", m=imputation_sum, lags="GPP")) ## lags by 1 day ##

  ## extract a nested list of data frames that just has the imputations###
  imputations_only <- map(full_amelia , ~.[["imputations"]])

  # Create empty lists to hold means and standard errors for parameter estimates across all the imputations    
  model_param_list_sim <- list()
  model_error_list_sim <- list()
  
  # Loop through each imputation and fit each with Arima, extracting estimates and errors
  for (i in seq_along(imputations_only)) {
    imputation_coefficients <- list()
    imputation_ses <- list()
    for (j in seq_along(imputations_only[[i]])) {
      imputation_fit <- Arima(imputations_only[[i]][[j]]$GPP, order = c(1,0,0), xreg = matrix(c(sim_missing_df [[j]][["light"]], sim_missing_df [[j]][["discharge"]]), ncol = 2))
      arima_coefs <- imputation_fit$coef
      arima_ses <- sqrt(diag(vcov(imputation_fit)))
      name <- paste('imp',seq_along((imputations_only)[[i]])[[j]],sep='')
      imputation_coefficients[[name]] <- arima_coefs
      imputation_ses[[name]] <- arima_ses
    }
    name1 <- names(imputations_only)[[i]]
    model_param_list_sim[[name1]] <- imputation_coefficients
    model_error_list_sim[[name1]] <- imputation_ses
  }
      
      
  ### Averages the models together back to 1 model per missing data prop ##
        
  coef_se_sim_list <- mapply(function(X,Y) {
    list(mi.meld(data.frame(X), data.frame(Y), byrow = FALSE))
  }, X = model_param_list_sim, Y = model_error_list_sim)
       
  return(list(paramlistsim = map(coef_se_sim_list , ~.["q.mi"]),
              selistsim = map(coef_se_sim_list , ~.["se.mi"])))
}


