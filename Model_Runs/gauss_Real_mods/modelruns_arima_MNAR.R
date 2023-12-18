# Load packages ## 
library(tidyverse)
#library(forecast) ## it hates this package...run with lowercase arima# 
library(Amelia)

## read in the data ##

gauss_real_MinMaxMiss <- readRDS("./data/missingDatasets/gauss_real_MinMaxMiss.rds")
pine_river_full <- read_csv("data/pine_river_data_prepped.csv")

  CurSim <- 1
   
  
  fit_arima_dropmissing <- function(sim_list, sim_pars, 
                                    forecast = TRUE, forecast_days = 31,
                                    dat_full){
    
    simmissingdf <-lapply(X = sim_list, 
                          FUN = function(X) cbind.data.frame(GPP = X, 
                                                             light = sim_pars$light, 
                                                             discharge = sim_pars$Q)) # Q is discharge
    ## holdout data for forecasting
    if(forecast){
      simmissingdf <- lapply(simmissingdf, function(df) {
        df[1:(366-forecast_days), ]  # Remove to save these for forecasting
      })
    }
    ## drop the missing values ###
    sim_missing_list_drop <- lapply(seq_along(simmissingdf), function(j) {
      drop_na(simmissingdf[[j]])
    })
    
    # fit arima models to list of datasets
    Arimaoutputdrop <- lapply(seq_along(sim_missing_list_drop ), function(j) {
      xreg1<-sim_missing_list_drop [[j]][["light"]]
      xreg2<-sim_missing_list_drop [[j]][["Q"]]
      modeldrop <- Arima(sim_missing_list_drop [[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(xreg1,xreg2), ncol = 2))
      arimacoefsdrop <-c(modeldrop$coef, modeldrop$sigma2)
      names(arimacoefsdrop) <- c("ar1", "intercept", "xreg1", "xreg2", "sigma")
      arimasesdrop<-sqrt(diag(vcov(modeldrop)))
      names(arimasesdrop) <- c("ar1", "intercept", "xreg1", "xreg2")
      arima_pars_j <- data.frame("parameters" = names(arimacoefsdrop), 
                                 "param_value" = arimacoefsdrop, 
                                 "param_se" = c(arimasesdrop,NA))
      outList <- list(arima_model = modeldrop,
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
        select(date, GPP, light, Q) %>% 
        rename("xreg1" = "light", "xreg2" = "Q")
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
  
  # Drop missing (complete case) + arima function ---------------------------------------------------
  
  fit_arima_dropmissing_CC <- function(sim_list, sim_pars, 
                                       forecast = TRUE, forecast_days = 31,
                                       dat_full){
    
    simmissingdf <-lapply(X = sim_list, 
                          FUN = function(X) cbind.data.frame(GPP = X, 
                                                             light = sim_pars$light, 
                                                             discharge = sim_pars$Q)) # Q is discharge
    
    ## holdout data for forecasting
    if(forecast){
      simmissingdf <- lapply(simmissingdf, function(df) {
        df[1:(366-forecast_days), ]  # Remove to save these for forecasting
      })
    }
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
      xreg2<-sim_missing_list_drop [[j]][["Q"]]
      modeldrop <- arima(sim_missing_list_drop [[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(xreg1,xreg2), ncol = 2))
      arimacoefsdrop <-c(modeldrop$coef, modeldrop$sigma2)
      names(arimacoefsdrop) <- c("ar1", "intercept", "xreg1", "xreg2", "sigma")
      arimasesdrop<-sqrt(diag(vcov(modeldrop)))
      names(arimasesdrop) <- c("ar1", "intercept", "xreg1", "xreg2")
      arima_pars_j <- data.frame("parameters" = names(arimacoefsdrop), 
                                 "param_value" = arimacoefsdrop, 
                                 "param_se" = c(arimasesdrop,NA))
      outList <- list(arima_model = modeldrop,
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
        select(date, GPP, light, Q) %>% 
        rename(xreg1 = "light", xreg2 = "Q")
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
  
  # arima + Kalman filter function ------------------------------------------
  
  fit_arima_Kalman <- function(sim_list, sim_pars, forecast = TRUE, forecast_days = 31,
                               dat_full){
    
    simmissingdf <-lapply(X = sim_list, 
                          FUN = function(X) cbind.data.frame(GPP = X, 
                                                             light = sim_pars$light, 
                                                             discharge = sim_pars$Q))
    ## holdout data for forecasting
    if(forecast){
      simmissingdf <- lapply(simmissingdf, function(df) {
        df[1:(366-forecast_days), ]  # Remove to save these for forecasting
      })
    }
    
    ## fit ARIMA with the missing values as NAS . Applies KALMAN FILTER###
    ArimaoutputNAs <- lapply(seq_along(simmissingdf), function(j) {
      xreg1<-simmissingdf [[j]][["light"]]
      xreg2<-simmissingdf [[j]][["Q"]]
      modelNAs <- arima(simmissingdf[[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(xreg1,xreg2), ncol = 2))
      arimacoefsNAs <- c(modelNAs$coef, modelNAs$sigma2)
      names(arimacoefsNAs) <- c("ar1", "intercept", "xreg1", "xreg2", "sigma")
      arimasesNAs<-sqrt(diag(vcov(modelNAs)))
      names(arimasesNAs) <- c("ar1", "intercept", "xreg1", "xreg2")
      arima_pars_j <- data.frame("parameters" = names(arimacoefsNAs), 
                                 "param_value" = arimacoefsNAs, 
                                 "param_se" = c(arimasesNAs,NA))
      outList <- list(arima_model = modelNAs,
                      arima_pars = arima_pars_j,
                      sim_params = sim_pars)
      
      return(outList)
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
        select(date, GPP, light, Q) %>% 
        rename(xreg1 = "light", xreg2 = "Q")
      xreg1 <- dat_forecast$xreg1
      xreg2 <- dat_forecast$xreg2
      predictions <- map_df(ArimaoutputNAs, function(mod){
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
  
  # multiple imputations + Arima function -----------------------------------
  
  fit_arima_MI <- function(sim_list, sim_pars, imputationsnum, forecast = TRUE, forecast_days = 31,
                           dat_full){
    
    days<-seq(1, 366)
    
    simmissingdf <-lapply(X = sim_list, 
                          FUN = function(X) cbind.data.frame(days= days,
                                                             GPP = X, 
                                                             light = sim_pars$light, 
                                                             discharge = sim_pars$Q))
    ## holdout data for forecasting
    if(forecast){
      simmissingdf <- lapply(simmissingdf, function(df) {
        df[1:(366-forecast_days), ]  # Remove to save these for forecasting
      })
    }
    
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
        xreg2<-amelias11sim [[i]][[j]][["Q"]]
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
    
    paramlistsim <- map_df(seq(1:16),
                           function(x){
                             data.frame("parameters" = c("intercept", "xreg1", "xreg2", "phi", "sigma"),
                                        "param_value" = c(listcoefsessim[[x]]$q.mi, sigmas[x]),
                                        "param_se" = c(listcoefsessim[[x]]$se.mi,NA)) 
                           }, 
                           .id = "tempNum")
    # update missingPropAutocor column
    numName_df <- data.frame("tempNum" = c(1:16), 
                             "missingprop_autocor" = names(listcoefsessim)) %>% 
      mutate("tempNum" = as.character(tempNum))
    
    paramlistsim <- left_join(paramlistsim, numName_df) %>% 
      select(missingprop_autocor, parameters, param_value, param_se)
    
    selistsim <- lapply(listcoefsessim, function(x) x$se.mi)
    
    # reframe list of coefficients and s.e.s into a single object for forecasting
    
    forecastList <- map(c(1:16), function(x) {
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
        select(date, GPP, light, Q) %>% 
        rename(xreg1 = "light", xreg2 = "Q")
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
  
  
  # Run models with drop missing + arima ------------------------------------
  
  arima_drop_MAR <- fit_arima_dropmissing(gauss_real_MinMaxMiss[[CurSim]]$y,gauss_real_MinMaxMiss[[CurSim]]$sim_params, forecast = TRUE, forecast_days = 31, dat_full = pine_river_full)
  
  ## formatting for figure
  # save arima model parameters
  arimadrop_MAR_df <- arima_drop_MAR$arima_pars
  arimadrop_MAR_df$missingness <- 'MAR'
  arimadrop_MAR_df$type <- 'dropNA_simple'
  arimadrop_MAR_df$run_no <- CurSim
  
  # save arima forecasts
  arimadrop_MAR_preds <- arima_drop_MAR$arima_forecast
  arimadrop_MAR_preds$missingness <- 'MAR'
  arimadrop_MAR_preds$type <- 'dropNA_simple'
  arimadrop_MAR_preds$run_no <- CurSim
  
  # Run models with drop missing complete case + arima ------------------------------------
  
  arima_dropCC_MAR <- fit_arima_dropmissing_CC(gauss_real_MinMaxMiss[[CurSim]]$y,gauss_real_MinMaxMiss[[CurSim]]$sim_params, forecast = TRUE, forecast_days = 31, dat_full = pine_river_full)
  
  ## formatting for figure
  # save arima model parameters
  arimadropCC_MAR_df <- arima_dropCC_MAR$arima_pars
  arimadropCC_MAR_df$missingness <- 'MAR'
  arimadropCC_MAR_df$type <- 'dropNA_complete'
  arimadropCC_MAR_df$run_no <- CurSim
  
  # save arima forecasts
  arimadropCC_MAR_preds <- arima_dropCC_MAR$arima_forecast
  arimadropCC_MAR_preds$missingness <- 'MAR'
  arimadropCC_MAR_preds$type <- 'dropNA_complete'
  arimadropCC_MAR_preds$run_no <- CurSim
  
  # Run models w/ Kalman filter + arima fxn ---------------------------------
  
  arima_kalman_MAR<- fit_arima_Kalman(gauss_real_MinMaxMiss[[CurSim]]$y,gauss_real_MinMaxMiss[[CurSim]]$sim_params, 
                                      forecast = TRUE, forecast_days = 31, dat_full = pine_river_full)
  ## formatting for figure
  # save arima model parameters
  arimaKalman_MAR_df <- arima_kalman_MAR$arima_pars
  arimaKalman_MAR_df$missingness <- 'MAR'
  arimaKalman_MAR_df$type <- 'Kalman Filter'
  arimaKalman_MAR_df$run_no <- CurSim
  
  # save arima forecasts
  arimaKalman_MAR_preds <- arima_kalman_MAR$arima_forecast
  arimaKalman_MAR_preds$missingness <- 'MAR'
  arimaKalman_MAR_preds$type <- 'Kalman Filter'
  arimaKalman_MAR_preds$run_no <- CurSim
  
  # Run models w/ Multiple Imputations --------------------------------------
  
  arima_mi_MAR <-  fit_arima_MI(gauss_real_MinMaxMiss[[CurSim]]$y,gauss_real_MinMaxMiss[[CurSim]]$sim_params, imputationsnum=5,
                                forecast = TRUE, forecast_days = 31,
                                dat_full = pine_river_full)
  ## formatting for figure
  # save arima model parameters
  arimaMI_MAR_df <- arima_mi_MAR$arima_pars
  arimaMI_MAR_df$missingness <- 'MAR'
  arimaMI_MAR_df$type <- 'Multiple Imputations'
  arimaMI_MAR_df$run_no <- CurSim
  
  # save arima forecasts
  arimaMI_MAR_preds <- map_df(arima_mi_MAR$arima_forecast, ~as.data.frame(.x),
                              .id = "missingprop_autocor")
  arimaMI_MAR_preds$missingness <- 'MAR'
  arimaMI_MAR_preds$type <- 'Multiple Imputations'
  arimaMI_MAR_preds$run_no <- CurSim
  
  # Combine output data and save --------------------------------------------
  
  paramarimaall_df<-rbind(arimadrop_MAR_df, arimadropCC_MAR_df, arimaKalman_MAR_df, arimaMI_MAR_df)
  
  Output <- matrix(data=NA, nrow=nrow(paramarimaall_df), ncol=ncol(paramarimaall_df))
  
  # Save the results of the current script's simulation to the appropriate column of output
  Output<- paramarimaall_df
  
  # add all the output data together
  valsAll <-cbind(CurSim,Output)
  
  # Write the output to the folder which will contain all output files as separate csv
  #    files with a single line of data.

  ## save prediction data
  
  paramarimaall_preds<-rbind(arimadrop_MAR_preds, arimadropCC_MAR_preds, arimaKalman_MAR_preds, arimaMI_MAR_preds)
  
  Output3 <- matrix(data=NA, nrow=nrow(paramarimaall_preds), ncol=ncol(paramarimaall_preds))
  
  # Save the results of the current script's simulation to the appropriate column of output
  Output3<- paramarimaall_preds
  
  # add all the output data together
  predsAll <-cbind(CurSim,Output3)
  
  # Write the output to the folder which will contain all output files as separate csv
  #    files with a single line of data.



write.csv(predsAll, file = "./data/model_results/gauss_real_MNAR_arima_FORECASTpreds.csv")

write.csv(valsAll, file = "./data/model_results/gauss_real_MNAR_arima_FORECASTvals.csv")