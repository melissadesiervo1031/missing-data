# Load packages
library(here)
library(tidyverse)
library(lubridate)

#' Function that will drop missing values and then fit model using ARIMA
#'
#' @param sim_list a list of GPP data from multiple simulations
#' @param sim_pars a list of simulation parameters corresponding to the GPP simulated values
#'
#' @return List of the ARIMA parameter estimates, errors, and true values from the simulations.
#'
#' @examples
#' 
#' gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss.rds")
#' fit_ricker_cc(y)
#' GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]
#' 
#' arima_drop <- fit_arima_dropmissing(GPP_sim_MAR$y,GPP_sim_MAR$sim_params)
#' 
fit_arima_dropmissing <- function(sim_list, sim_pars){
  
  simmissingdf <-lapply(X = sim_list, 
                        FUN = function(X) cbind.data.frame(GPP = X, 
                                                           light = sim_pars$X[,2], 
                                                           discharge = sim_pars$X[,3]))
  ## drop the missing values ###
  
  sim_missing_list_drop <- lapply(seq_along(simmissingdf), function(j) {
    drop_na(simmissingdf[[j]])
  })
  
  # fit arima models to list of datasets
  Arimaoutputdrop <- lapply(seq_along(sim_missing_list_drop ), function(j) {
    modeldrop <- Arima(sim_missing_list_drop [[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(sim_missing_list_drop [[j]][["light"]],sim_missing_list_drop [[j]][["discharge"]]), ncol = 2))
    arimacoefsdrop<-modeldrop$coef
    arimasesdrop<-sqrt(diag(vcov(modeldrop)))
  
    return(list(arima_pars = arimacoefsdrop,
                arima_errors = arimasesdrop,
                sim_params = sim_pars))
        })
}

