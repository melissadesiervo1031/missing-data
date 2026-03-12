# Load packages ## 
#make sure these are already in the folder on supercomputer where I need them ##

#.libPaths("/pfs/tc1/home/astears/R/x86_64-pc-linux-gnu-library/4.2")

library(tidyverse)
#library(forecast) ## it hates this package...run with lowercase arima# 
library(Amelia)

# This script will run 3 ARIMA functions (drop missing, Kalman, Multiple imputations 
#over a nested list with increasing prop missing, over 1000+ simulations ###)

#CurSim = like a loop ##

# CurSim <- commandArgs(trailingOnly = TRUE) #Look at command line arguments only after the R script
# CurSim <- as.numeric(CurSim)
# CurSim <- CurSim + 1 # since the Slurm array is 0 indexed

## read in the autocor_01 list ##

#gauss_sim_randMiss_autoCorr_01 <- readRDS("/project/modelscape/users/astears/gauss_sim_randMiss_B.rds")
gauss_sim_randMiss_autoCorr_01 <- readRDS("./data/missingDatasets/gauss_sim_randMiss_B.rds")
# make file for output beforehand in supercomputer folder 
# will put them all together after all run, using the command line
#OutFile <- paste("gauss_sim_randMiss_modResults_A/", CurSim, "arimavals.csv", sep = "")
OutFile <- paste("./data/model_results/gauss_sim_randMiss_modelResults_B/")
#########################################################################################
### MY ARIMA FUNCTIONS #####
##########################################################################################

fit_arima_dropmissing <- function(sim_list, sim_pars, 
                                  forecast = TRUE, forecast_days = 73,
                                  dat_full){
  
  simmissingdf <-lapply(X = sim_list, 
                        FUN = function(X) cbind.data.frame(GPP = X[1:292], #the first element of the list, which includes the full length of the time series w/ no msisingness, needs to be curtailed to match the lenght of the other missing datasets 
                                                           light = sim_pars$X[,2][1:292], 
                                                           discharge = sim_pars$X[,3][1:292])) # Q is discharge
  ## data for forecasting has arleady been held out 
  
  ## drop the missing values ###
  sim_missing_list_drop <- lapply(seq_along(simmissingdf), function(j) {
    drop_na(simmissingdf[[j]])
  })
  
  # fit arima models to list of datasets
  Arimaoutputdrop <- lapply(seq_along(sim_missing_list_drop ), function(j) {
    xreg1<-sim_missing_list_drop [[j]][["light"]]
    xreg2<-sim_missing_list_drop [[j]][["discharge"]]
    modeldrop <- arima(sim_missing_list_drop [[j]][["GPP"]], order = c(1,0,0), xreg = matrix(c(xreg1,xreg2), ncol = 2))
    arimacoefsdrop <-c(modeldrop$coef, modeldrop$sigma2)
    names(arimacoefsdrop) <- c("ar1", "intercept", "xreg1", "xreg2", "sigma")
    arimasesdrop<-sqrt(diag(vcov(modeldrop)))
    names(arimasesdrop) <- c("ar1", "intercept", "xreg1", "xreg2")
    arima_pars_j <- data.frame("parameters" = names(arimacoefsdrop), 
                               "param_value" = arimacoefsdrop, 
                               "param_se" = c(arimasesdrop,NA))
    outList <- list(arima_model = modeldrop,
                    arima_errors = arimasesdrop,
                    arima_pars = arima_pars_j,
                    sim_params = sim_pars)
    
    return(outList)
  })
  names(Arimaoutputdrop) <- names(simmissingdf)
  # flatten arima parameters
  arima_pars <- map_df(Arimaoutputdrop,
                       function(x) {
                         x[["arima_pars"]]
                       }, 
                       .id = "missingprop_autocor")
  rownames(arima_pars) <- NULL
  
  if(forecast){  
    dat_forecast <- dat_full %>%
      slice((nrow(dat_full)-forecast_days):nrow(dat_full)) %>%
      select(date, GPP, light, discharge) %>% 
      rename("xreg1" = "light", "xreg2" = "discharge")
    xreg1 <- dat_forecast$xreg1
    xreg2 <- dat_forecast$xreg2
    predictions <- map_df(Arimaoutputdrop, function(mod){
      data.frame(predict(mod$arima_model, n.ahead = forecast_days+1, 
                         newxreg = dat_forecast[,c(3,4)]),
                 "date" = dat_forecast$date,
                 "GPP" = dat_forecast$GPP)
    },.id = "missingprop_autocor") 
    
    return(list(arima_forecast = predictions,
                arima_pars = arima_pars,
                sim_params = sim_pars))
  }
  
}
# fit complete case drop missing

fit_arima_dropmissing_CC <- function(sim_list, sim_pars, 
                                     forecast = TRUE, forecast_days = 73,
                                     dat_full){
  
  simmissingdf <- lapply(X = sim_list, 
                         FUN = function(X) cbind.data.frame(GPP = X[1:292], #the first element of the list, which includes the full length of the time series w/ no msisingness, needs to be curtailed to match the lenght of the other missing datasets 
                                                            light = sim_pars$X[,2][1:292], 
                                                            discharge = sim_pars$X[,3][1:292])) # Q is discharge
  
  ## data have already been held out for forecasting
  
  # remove data in a "complete case" way
  # compile into sliced dataframe
  sim_missing_list_drop <- map(simmissingdf, function(x) {
    temp <- data.frame(
      yt = x[2:nrow(x),],
      ytm1 = x[1:(nrow(x)-1),]
    )
    # drop incomplete cases
    x[complete.cases(temp),]
  }
  )
  
  
  # fit arima models to list of datasets
  Arimaoutputdrop <- lapply(seq_along(sim_missing_list_drop ), function(j) {
    xreg1<-sim_missing_list_drop [[j]][["light"]]
    xreg2<-sim_missing_list_drop [[j]][["discharge"]]
    modeldrop <- arima(sim_missing_list_drop [[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(xreg1,xreg2), ncol = 2))
    arimacoefsdrop <-c(modeldrop$coef, modeldrop$sigma2)
    names(arimacoefsdrop) <- c("ar1", "intercept", "xreg1", "xreg2", "sigma")
    arimasesdrop<-sqrt(diag(vcov(modeldrop)))
    names(arimasesdrop) <- c("ar1", "intercept", "xreg1", "xreg2")
    arima_pars_j <- data.frame("parameters" = names(arimacoefsdrop), 
                               "param_value" = arimacoefsdrop, 
                               "param_se" = c(arimasesdrop,NA))
    outList <- list(arima_model = modeldrop,
                    arima_errors = arimasesdrop,
                    arima_pars = arima_pars_j,
                    sim_params = sim_pars)
    
    return(outList)
  })
  
  names(Arimaoutputdrop) <- names(simmissingdf)
  
  # flatten arima parameters
  arima_pars <- map_df(Arimaoutputdrop,
                       function(x) {
                         x[["arima_pars"]]
                       }, 
                       .id = "missingprop_autocor")
  rownames(arima_pars) <- NULL
  
  if(forecast){  
    dat_forecast <- dat_full %>%
      slice((nrow(dat_full)-forecast_days):nrow(dat_full)) %>%
      select(date, GPP, light, discharge) %>% 
      rename(xreg1 = "light", xreg2 = "discharge")
    xreg1 <- dat_forecast$xreg1
    xreg2 <- dat_forecast$xreg2
    predictions <- map_df(Arimaoutputdrop, function(mod){
      data.frame(predict(mod$arima_model, n.ahead = forecast_days+1, 
                         newxreg = dat_forecast[,c(3,4)]),
                 "date" = dat_forecast$date,
                 "GPP" = dat_forecast$GPP)
    },.id = "missingprop_autocor") 
    
    return(list(arima_forecast = predictions,
                arima_pars = arima_pars,
                sim_params = sim_pars))
  }
  
}
### Function that will have missing values as NA and then fit model using ARIMA w/ Kalman filter ###

fit_arima_Kalman <- function(sim_list, sim_pars, forecast = TRUE, forecast_days = 73,
                             dat_full){
  
  simmissingdf <- lapply(X = sim_list, 
                         FUN = function(X) cbind.data.frame(GPP = X[1:292], #the first element of the list, which includes the full length of the time series w/ no msisingness, needs to be curtailed to match the lenght of the other missing datasets 
                                                            light = sim_pars$X[,2][1:292], 
                                                            discharge = sim_pars$X[,3][1:292])) # Q is discharge
  
  
  ## fit ARIMA with the missing values as NAS . Applies KALMAN FILTER###
  ArimaoutputNAs <- lapply(seq_along(simmissingdf), function(j) {
    xreg1<-simmissingdf [[j]][["light"]]
    xreg2<-simmissingdf [[j]][["discharge"]]
    modelNAs <- try(arima(simmissingdf[[j]][["GPP"]], order = c(1,0,0), xreg = matrix(c(xreg1,xreg2), ncol = 2)))
    if (is.list(modelNAs)) {
      arimacoefsNAs <- c(modelNAs$coef, modelNAs$sigma2)
      names(arimacoefsNAs) <- c("ar1", "intercept", "xreg1", "xreg2", "sigma")
      arimasesNAs<-sqrt(diag(vcov(modelNAs)))
      names(arimasesNAs) <- c("ar1", "intercept", "xreg1", "xreg2")
      arima_pars_j <- data.frame("parameters" = names(arimacoefsNAs), 
                                 "param_value" = arimacoefsNAs, 
                                 "param_se" = c(arimasesNAs,NA))
      outList <- list(arima_model = modelNAs,
                      arima_errors = arima_pars_j$param_se,
                      arima_pars = arima_pars_j,
                      sim_params = sim_pars)
      return(outList)
    } else {
      arimacoefsNAs <- c(NA, NA, NA, NA, NA)
      names(arimacoefsNAs) <- c("ar1", "intercept", "xreg1", "xreg2", "sigma")
      arimasesNAs<- c(NA, NA, NA, NA)
      names(arimasesNAs) <- c("ar1", "intercept", "xreg1", "xreg2")
      arima_pars_j <- data.frame("parameters" = names(arimacoefsNAs), 
                                 "param_value" = arimacoefsNAs, 
                                 "param_se" = c(arimasesNAs,NA))
      outList <- list(arima_model = "modelFitWasSingular",
                      arima_errors = arima_pars_j$param_se,
                      arima_pars = arima_pars_j,
                      sim_params = sim_pars)
      return(outList)
    }
    
  })
  
  
  names(ArimaoutputNAs) <- names(simmissingdf)
  
  # flatten arima parameters
  arima_pars <- map_df(ArimaoutputNAs,
                       function(x) {
                         x[["arima_pars"]]
                       }, 
                       .id = "missingprop_autocor")
  rownames(arima_pars) <- NULL
  
  if(forecast){  
    dat_forecast <- dat_full %>%
      slice((nrow(dat_full)-forecast_days):nrow(dat_full)) %>%
      select(date, GPP, light, discharge) %>% 
      rename(xreg1 = "light", xreg2 = "discharge")
    xreg1 <- dat_forecast$xreg1
    xreg2 <- dat_forecast$xreg2
    predictions <- map_df(ArimaoutputNAs[1], function(mod){
      if (!is.list(mod$arima_model)) {
        data.frame( "pred" = rep.int(NA, times = forecast_days+1),
                    "se" = rep.int(NA, times = forecast_days+1),
                    "date" = rep.int(NA, times = forecast_days+1), "GPP" = rep.int(NA, times = forecast_days+1))
      } else {
        data.frame(predict(mod$arima_model, n.ahead = forecast_days+1, 
                           newxreg = dat_forecast[,c(3,4)]),
                   "date" = dat_forecast$date,
                   "GPP" = dat_forecast$GPP)
      }
      
    },.id = "missingprop_autocor") 
    
    return(list(arima_forecast = predictions,
                arima_pars = arima_pars,
                sim_params = sim_pars))
  }
}

###### 

### Function that will impute missing values w/ AMELIA and then fit model using ARIMA ###


fit_arima_MI <- function(sim_list, sim_pars, imputationsnum, forecast = TRUE, forecast_days = 73,
                         dat_full){
  
  simmissingdf <- lapply(X = sim_list, 
                         FUN = function(X) cbind.data.frame(days = 1:292,
                                                            GPP = X[1:292], #the first element of the list, which includes the full length of the time series w/ no msisingness, needs to be curtailed to match the lenght of the other missing datasets 
                                                            light = sim_pars$X[,2][1:292], 
                                                            discharge = sim_pars$X[,3][1:292])) # Q is discharge
  
  
  amelia1sim <-lapply(X = simmissingdf  , FUN = function(X)   amelia(X, ts="days", m=imputationsnum, lags="GPP")) ## lags by 1 day ##
  
  
  ##nested list of dataframes that just has the imputations###
  amelias11sim<-map(amelia1sim , ~.[["imputations"]])
  
  ##forloop## gives us the model parameters and errors for all the ARIMA models on imputed datasets
  
  modelparamlistsim <- list()
  modelerrorlistsim <- list()
  modelobjectlist <- list()
  for (i in seq_along(amelias11sim)) {
    a=list()
    aa=list()
    mod_a <- list()
    for (j in seq_along(amelias11sim[[i]])) {
      xreg1<-amelias11sim [[i]][[j]][["light"]]
      xreg2<-amelias11sim [[i]][[j]][["discharge"]]
      tempobj=arima(amelias11sim[[i]][[j]]$GPP, order = c(1,0,0), xreg = matrix(c(xreg1, xreg2), ncol = 2))
      arimacoefs<-c(tempobj$coef, tempobj$sigma2)
      names(arimacoefs) <- c("ar1", "intercept", "xreg1", "xreg2", "sigma")
      arimases<-sqrt(diag(vcov(tempobj)))
      names(arimases) <- c("ar1", "intercept", "xreg1", "xreg2")
      name <- paste('imp',seq_along((amelias11sim)[[i]])[[j]],sep='')
      a[[name]] <- arimacoefs
      aa[[name]]<-arimases
      mod_a[[name]] <- tempobj
    }
    #name1 <- names(amelias11sim)[[i]]
    modelparamlistsim[[i]] <- a
    modelerrorlistsim[[i]] <- aa
    modelobjectlist[[i]] <- mod_a
  }
  
  ### Averages the models together back to 1 model per missing data prop ##
  
  listcoefsessim<-mapply(function(X,Y) {
    list(mi.meld(data.frame(X), data.frame(Y), byrow=FALSE))
  }, X=lapply(modelparamlistsim, function(x) lapply(x, function (x) x[1:4])), 
  Y=modelerrorlistsim)
  # rename coef list
  names(listcoefsessim) <- names(amelias11sim)
  
  # put the mean sigma for each imputation back into the listcoefsessim list
  sigmas_temp <- sapply(modelparamlistsim, function(x) lapply(x, function(x) x[5]))
  sigmas <- apply(sigmas_temp, MARGIN = 2, function(x) mean(as.numeric(x), na.rm = TRUE))
  
  # make return values
  #paramlistsim <- map(listcoefsessim , ~.["q.mi"])
  
  paramlistsim <- map_df(seq(1:length(listcoefsessim)),
                         function(x){
                           data.frame("parameters" = c("phi", "intercept", "xreg1", "xreg2",  "sigma"),
                                      "param_value" = c(listcoefsessim[[x]]$q.mi, sigmas[x]),
                                      "param_se" = c(listcoefsessim[[x]]$se.mi,NA)) 
                         }, 
                         .id = "tempNum")
  # update missingPropAutocor column
  numName_df <- data.frame("tempNum" = c(1:length(listcoefsessim)), 
                           "missingprop_autocor" = names(listcoefsessim)) %>% 
    mutate("tempNum" = as.character(tempNum))
  
  paramlistsim <- left_join(paramlistsim, numName_df) %>% 
    select(missingprop_autocor, parameters, param_value, param_se)
  
  selistsim <- lapply(listcoefsessim, function(x) x$se.mi)
  
  # reframe list of coefficients and s.e.s into a single object for forecasting
  
  forecastList <- map(c(1:length(listcoefsessim)), function(x) {
    test <- modelobjectlist[[i]]$imp1  
    test$coef <- as.vector(listcoefsessim[[x]]$q.mi)
    names(test$coef) <- c("ar1","intercept","matrix(c(xreg1, xreg2), ncol = 2)1","matrix(c(xreg1, xreg2), ncol = 2)2")
    test$sigma2 <- sigmas[[i]]
    return(test)
  }
  ) 
  names(forecastList) <- names(listcoefsessim) 
  
  if(forecast){  
    dat_forecast <- dat_full %>%
      slice((nrow(dat_full)-forecast_days):nrow(dat_full)) %>%
      select(date, GPP, light, discharge) %>% 
      rename(xreg1 = "light", xreg2 = "discharge")
    xreg1 <- dat_forecast$xreg1
    xreg2 <- dat_forecast$xreg2
    predictions <- map_df(forecastList, function(mod){
      data.frame(predict(mod, n.ahead = forecast_days+1, 
                         newxreg = dat_forecast[,c(3,4)]),
                 "date" = dat_forecast$date,
                 "GPP" = dat_forecast$GPP)
    },.id = "missingprop_autocor") 
    
    
    return(list(arima_forecast = predictions,
                arima_pars = paramlistsim,
                arima_se = selistsim,
                sim_params = sim_pars))
  }
  
  return(list(paramlistsim, 
              selistsim
  ))
}


for (i in 1:5000) {
  CurSim <- i
  #### set up data for this iteration
  
  gauss_sim_CurSim_vals <- gauss_sim_randMiss_autoCorr_01[[CurSim]]$y
  #gauss_sim_CurSim_vals <- gauss_sim_randMiss_autoCorr_01$gauss_sim997_randMiss_autoCorr_50$y
  gauss_sim_CurSim_params <- gauss_sim_randMiss_autoCorr_01[[CurSim]]$sim_params
  #gauss_sim_CurSim_params <- gauss_sim_randMiss_autoCorr_01$gauss_sim997_randMiss_autoCorr_50$sim_params
  
  # "Full" dataset (i.e. the first element in the gauss_sim_CurSim_vals list)
  dat_full <- data.frame("date" = 1:365,
                         "GPP" = gauss_sim_CurSim_vals$y_noMiss,
                         "light" = gauss_sim_CurSim_params$X[,2],
                         "discharge" = gauss_sim_CurSim_params$X[,3]
  )
  #####################################################
  #### MODEL RUN ARIMA DROP ##############
  #########################################################
  
  ## for this model-type only, fit a model on the full dataset and forecast with 
  # it (for all others, only fit models and predict w/ missing datasets)
  arima_drop_MAR<- fit_arima_dropmissing(sim_list = gauss_sim_CurSim_vals,
                                         sim_pars = gauss_sim_CurSim_params, 
                                         dat_full = dat_full, forecast = TRUE, forecast_days = 73)
  
  
  ########### formatting for figure #############
  ## formatting for figure
  # save arima model parameters
  arimadrop_MAR_df <- arima_drop_MAR$arima_pars
  arimadrop_MAR_df$missingness <- 'MAR'
  arimadrop_MAR_df$type <- 'dropNA_simple'
  arimadrop_MAR_df$curSim <- CurSim
  
  # save arima forecasts
  arimadrop_MAR_preds <- arima_drop_MAR$arima_forecast
  arimadrop_MAR_preds$missingness <- 'MAR'
  arimadrop_MAR_preds$type <- 'dropNA_simple'
  arimadrop_MAR_preds$curSim <- CurSim
  
  
  #####################################################
  #### MODEL RUN ARIMA DROP --complete case ##############
  #########################################################
  
  arima_drop_CC_MAR <- fit_arima_dropmissing_CC(sim_list = gauss_sim_CurSim_vals[2:16],
                                                sim_pars = gauss_sim_CurSim_params, 
                                                forecast = TRUE, 
                                                forecast_days = 73,
                                                dat_full = dat_full)
  
  
  ########### formatting for figure #############
  # save arima model parameters
  arimadropCC_MAR_df <- arima_drop_CC_MAR$arima_pars
  arimadropCC_MAR_df$missingness <- 'MAR'
  arimadropCC_MAR_df$type <- 'dropNA_complete'
  arimadropCC_MAR_df$curSim <- CurSim
  
  # save arima forecasts
  arimadropCC_MAR_preds <- arima_drop_CC_MAR$arima_forecast
  arimadropCC_MAR_preds$missingness <- 'MAR'
  arimadropCC_MAR_preds$type <- 'dropNA_complete'
  arimadropCC_MAR_preds$curSim <- CurSim
  
  #####################################################
  #### MODEL RUN KALMAN FILTER ##############
  #########################################################
  arima_kalman_MAR<- fit_arima_Kalman(sim_list = gauss_sim_CurSim_vals[2:16],
                                      sim_pars = gauss_sim_CurSim_params, 
                                      forecast = TRUE, 
                                      forecast_days = 73,
                                      dat_full = dat_full)
  
  
  ## pull out and label what we need ###
  arimaKalman_MAR_df <- arima_kalman_MAR$arima_pars
  arimaKalman_MAR_df$missingness <- 'MAR'
  arimaKalman_MAR_df$type <- 'Kalman Filter'
  arimaKalman_MAR_df$curSim <- CurSim
  
  # save arima forecasts
  arimaKalman_MAR_preds <- arima_kalman_MAR$arima_forecast
  arimaKalman_MAR_preds$missingness <- 'MAR'
  arimaKalman_MAR_preds$type <- 'Kalman Filter'
  arimaKalman_MAR_preds$curSim <- CurSim
  #######################################################################################################
  
  
  #####################################################
  #### MODEL RUN MULTIPLE IMPUTATIONS  ##############
  #########################################################
  
  arima_mi_MAR <-  fit_arima_MI(sim_list = gauss_sim_CurSim_vals[2:16],
                                sim_pars = gauss_sim_CurSim_params, 
                                forecast = TRUE, 
                                forecast_days = 73,
                                imputationsnum=5,
                                dat_full = dat_full)
  
  ##pulls out parameters and ses ##
  # save arima model parameters
  arimaMI_MAR_df <- arima_mi_MAR$arima_pars
  arimaMI_MAR_df$missingness <- 'MAR'
  arimaMI_MAR_df$type <- 'Multiple Imputations'
  arimaMI_MAR_df$curSim <- CurSim
  
  # save arima forecasts
  arimaMI_MAR_preds <- arima_mi_MAR$arima_forecast
  arimaMI_MAR_preds$missingness <- 'MAR'
  arimaMI_MAR_preds$type <- 'Multiple Imputations'
  arimaMI_MAR_preds$curSim <- CurSim
  #############################################################################################################
  
  ###################################################
  #### COMBINE ALL MODEL RUNS AND SAVE #########
  #################################################
  
  ###
  # put all model parameters together
  params_MAR_all <- rbind(
    arimadrop_MAR_df, 
    arimadropCC_MAR_df, 
    arimaKalman_MAR_df, 
    arimaMI_MAR_df
  )
  # put all forecast predictions together
  preds_MAR_all <- rbind(
    arimadrop_MAR_preds, 
    arimadropCC_MAR_preds,
    arimaKalman_MAR_preds, 
    arimaMI_MAR_preds
  )
  
  # add in the "simulation number" for this iteration (which is stored in the name of the data list element)
  simName <- str_sub(string = names(gauss_sim_randMiss_autoCorr_01[CurSim]), 
                     start = str_locate_all(
                       pattern = "_", string  = names(gauss_sim_randMiss_autoCorr_01[CurSim])
                     )[[1]][1,1] + 4,
                     end = str_locate_all(
                       pattern = "_", string  = names(gauss_sim_randMiss_autoCorr_01[CurSim])
                     )[[1]][2,1]-1
  )
  
  params_MAR_all$sim_no <- simName
  preds_MAR_all$sim_no <- simName
  
  # Write the output to the folder which will contain all output files as separate csv
  #    files with a single line of data.
  write.csv(params_MAR_all, file = paste0(OutFile, CurSim,"_params.csv"), row.names = FALSE)
  write.csv(preds_MAR_all, file = paste0(OutFile, CurSim,"_predValues.csv"), row.names = FALSE)
  
}

#to aggregate predictions into one doc.
#     awk '(NR == 1) || (FNR > 1)' *_params.csv > AllParams_arima.csv
#     awk '(NR == 1) || (FNR > 1)' *_predValues.csv > AllPreds_arima.csv
