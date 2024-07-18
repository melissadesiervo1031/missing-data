# Load packages
library(here)
library(tidyverse)
library(lubridate)
library(forecast)

#' Function that will fit the model using ARIMA while dealing with missing 
#' values using the Kalman filter ARIMA
#'
#' @param sim_list a list of GPP data from multiple simulations
#' @param sim_pars a list of simulation parameters corresponding to the GPP simulated values
#'
#' @return List of the ARIMA parameter estimates, errors, and true values from the simulations.
#'
#' @examples
#' 
#' gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss.rds")
#' GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]
#' 
#' arima_kalman <- fit_arima_Kalman(GPP_sim_MAR$y,GPP_sim_MAR$sim_params)
#' 

### Function that will drop missing values and then fit model using ARIMA ###

fit_arima_Kalman <- function(sim_list, sim_pars){
  
  simmissingdf <-lapply(X = sim_list, 
                        FUN = function(X) cbind.data.frame(GPP = X, 
                                                           light = sim_pars$X[,2], 
                                                           discharge = sim_pars$X[,3]))
  
  ## fit ARIMA with the missing values as NAS . Applies KALMAN FILTER###
  
  ArimaoutputNAs <- lapply(seq_along(simmissingdf), function(j) {
    modelNAs <- Arima(simmissingdf[[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(simmissingdf [[j]][["light"]],simmissingdf [[j]][["discharge"]]), ncol = 2))
    arimacoefsNAs<-modelNAs$coef
    arimasesNAs<-sqrt(diag(vcov(modelNAs)))
    list(arimacoefsNAs=arimacoefsNAs, arimasesNAs=arimasesNAs)
    
    return(list(arima_pars = arimacoefsNAs,
                arima_errors = arimasesNAs,
                sim_params = sim_pars))
  })
}
  
