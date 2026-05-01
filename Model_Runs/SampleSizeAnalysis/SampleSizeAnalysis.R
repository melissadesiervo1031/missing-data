#//////////////////
# Sample Size Analysis: Are the decreases in parameter recovery with increased 
# missingness due to more missing data or smaller sample size?
# 14 April 2026
#//////////////////

# additional analysis to likely put in a supplement: fit models to simulated
#datasets that have the same number of observations total, but different
#proportions of missing values (i.e. datasets that are 100 obs. long w/ 10\%
#missingness (90 values) and datasets that are 113 observations long w/ 20\%
#missingness (~90 values)) to address the question of whether it's actually the
#amount of missing data or just smaller datasets that is driving patterns we
#observe in parameter recovery


# load packages ------------------------------------------------------------

library(tidyverse)
library(Amelia)
library(brms)

# load internal functions -------------------------------------------------

source("./Functions/missing_data_functions.R")

# simulate data -----------------------------------------------------------
# argument to run simulations or use previously-run simulation data
simulateOpt <- FALSE
# argument to run models 
modelOpt <- TRUE
## Simulate datasets of different length
# 100 datasets each of length 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200

# use code from ./Simulations/guassian_ar1_data_sims.R
#' n Integer number of observations
#' p Integer number of covariates (columns)
#' sigma2 Positive-valued scalar or vector of variances for the covariates
#' mu Scalar or vector of means for the columns
#'
rand_mod_matrix <- function(n, p, sigma2 = 1, mu = 0, intercept = TRUE){
  M <- matrix(
    data = runif(p^2, min = -1), nrow = p, ncol = p
  )
  # create orthogonal matrix
  L <- qr.Q(qr(M))
  
  # create diagonal matrix of sds or scales
  D <- diag(sigma2, nrow = p, ncol = p)
  
  # create covariance matrix
  Sigma <- t(L) %*% D %*% L
  if(length(mu) == 1){
    mu <- rep(mu, p)
  }
  
  X <- mvtnorm::rmvnorm(
    n = n,
    mean = mu,
    sigma = Sigma
  )
  
  if(intercept){
    return(cbind(rep(1, n), X))
  } else{
    return(X)
  }
  
}

##### simulate datasets
if (simulateOpt) {
  set.seed(235)
  # vector of simulation lengths
  n_options <- c(100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200)
  sim_list <- lapply(n_options, function(x) {
    # global parameters
    nsims <- 100 # number of simulations
    n <- x # number of simulated observations
    p <- 2 # number of covariates
    sde <- 1
    
    # parameters for each simulation
    params <- purrr::map(
      1:nsims,
      ~ list(
        phi = runif(1, 0, 0.8), ###phi gets weird too close to 1 ###
        beta = rnorm(p + 1),
        X = rand_mod_matrix(n = n, p = p)
      )
    )
    
    # simulated datasets
    sims <- purrr::map(
      params,
      ~ list(
        y = {as.double(arima.sim(
          n = n,
          model = list(ar = .x$phi),
          sd = sde
        )) + as.double(.x$X %*% .x$beta)},
        sim_params = .x
      )
    )
    # assign simulation names to each list element
    names(sims) <- paste0("N_",x,"_sim_",1:100)
    return(sims)
  })
  names(sim_list) <- paste0("N_",c(100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200))
  
  
  # save data
  saveRDS(sim_list,
          file = here::here("data/gauss_ar1_0miss_datasets_forSampleSizeSupplement.rds")
  )
  # sim_list <- readRDS("data/gauss_ar1_0miss_datasets_forSampleSizeSupplement.rds")
  
  # Introduce missingness into simulated data -------------------------------
  sim_list_miss <- lapply(1:11, function(y) {
    list_y <- sim_list[[y]]
    # now apply across all iterations of simulations for the current time series length
    temp_y <- lapply(list_y, function(z) {
      # calculate the appropriate % missingness for this iteration (which is the amount of missingness needed to get to 100 observations)
      propMiss_temp <- 1 - (100/(names(sim_list)[y] %>% str_extract(pattern = "[0-9]+") %>% as.numeric) )
      # create missing dataset for this time series
      if (propMiss_temp == 0) {
        list_z <- c(z,
                    list("propMissAct_0" = z$y)
        )
      } else {
        
        n <- length(z$y)
        while(n != 100) {
          missing_y <- makeMissing(timeSeries = z$y, 
                                   typeMissing = "random", 
                                   autoCorr = 0, # use fixed autocorrelation in this simulation
                                   propMiss = propMiss_temp
          )
          n <- sum(!is.na(missing_y[[1]]))
        }
        
        # save
        list_z <- c(z,
                    missing_y    
        )
      }
      
      list_z$missList_N <- sum(!is.na(list_z[[3]]))
      return(list_z)
    })
    return(temp_y)
  })
  #names(sim_list_miss) <- paste0("N_",c(100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200))
  sim_list_miss_short <- sim_list_miss %>% 
    list_flatten()
  
  # save data
  saveRDS(sim_list_miss_short,
          file = here::here("data/gauss_ar1_0miss_datasets_WithMissingness_forSampleSizeSupplement.rds")
  )
} else {
  sim_list <- readRDS(here::here("data/gauss_ar1_0miss_datasets_forSampleSizeSupplement.rds"))
  sim_list_miss_short <- readRDS(here::here("data/gauss_ar1_0miss_datasets_WithMissingness_forSampleSizeSupplement.rds"))
}

# Fit models to missing datasets ------------------------------------------
# start w/ just using kalman filter missing data method
# function for kalman filter
fit_arima_Kalman <- function(sim_list, sim_pars, forecast = TRUE, forecast_days = 73,
                             dat_full){
  
  simmissingdf <- list(cbind.data.frame(GPP = sim_list, 
                                        light = sim_pars$X[,2],
                                        discharge = sim_pars$X[,3]))
  # lapply(X = sim_list, 
  #                        FUN = function(X) cbind.data.frame(GPP = X[1:292], #the first element of the list, which includes the full length of the time series w/ no msisingness, needs to be curtailed to match the lenght of the other missing datasets 
  #                                                           light = sim_pars$X[,2][1:292], 
  #                                                           discharge = sim_pars$X[,3][1:292])) # Q is discharge
  
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
    predictions <- map_df(ArimaoutputNAs, function(mod){
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
  return(list(#arima_forecast = predictions,
    arima_pars = arima_pars,
    sim_params = sim_pars))
}

# fit simple case drop missing
fit_arima_dropmissing <- function(sim_list, sim_pars, 
                                  forecast = TRUE, forecast_days = 73,
                                  dat_full){
  
  simmissingdf <- list(cbind.data.frame(GPP = sim_list, 
                                        light = sim_pars$X[,2],
                                        discharge = sim_pars$X[,3]))
  ## data for forecasting has arleady been held out 
  
  ## drop the missing values ###
  sim_missing_list_drop <- lapply(seq_along(simmissingdf), function(j) {
    drop_na(simmissingdf[[j]])
  })
  
  # fit arima models to list of datasets
  Arimaoutputdrop <- lapply(seq_along(sim_missing_list_drop ), function(j) {
    xreg1<-sim_missing_list_drop [[j]][["light"]]
    xreg2<-sim_missing_list_drop [[j]][["discharge"]]
    modeldrop <- arima(sim_missing_list_drop [[j]][["GPP"]], order = c(1,0,0), xreg = matrix(c(xreg1,xreg2), ncol = 2) )
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
  
  return(list(#arima_forecast = predictions,
    arima_pars = arima_pars,
    sim_params = sim_pars))
  
}
# fit complete case drop missing
fit_arima_dropmissing_CC <- function(sim_list, sim_pars, 
                                     forecast = TRUE, forecast_days = 73,
                                     dat_full){
  
  simmissingdf <- list(cbind.data.frame(GPP = sim_list, 
                                        light = sim_pars$X[,2],
                                        discharge = sim_pars$X[,3]))  
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
    modeldrop <- arima(sim_missing_list_drop [[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(xreg1,xreg2), ncol = 2) )
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
  
  return(list(#arima_forecast = predictions,
    arima_pars = arima_pars,
    sim_params = sim_pars)
  )
  
}

# multiple imputation
fit_arima_MI <- function(sim_list, sim_pars, imputationsnum, forecast = TRUE, forecast_days = 73,
                         dat_full){
  
  simmissingdf <- list(cbind.data.frame(days = 1:length(sim_list),
                                        GPP = sim_list, 
                                        light = sim_pars$X[,2],
                                        discharge = sim_pars$X[,3]))  
  
  
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
                           "missingprop_autocor" =NA) %>% 
    mutate("tempNum" = as.character(tempNum))
  
  paramlistsim <- left_join(paramlistsim, numName_df) %>% 
    select(missingprop_autocor, parameters, param_value, param_se)
  
  selistsim <- lapply(listcoefsessim, function(x) x$se.mi)
  
  return(list(paramlistsim, 
              selistsim
  ))
}

# brms
fit_brms_model <- function(sim_list, sim_pars, 
                           iter = 4000, include_missing = FALSE,
                           forecast = TRUE, forecast_days = 73,
                           dat_full ){
  
  simmissingdf <- list(cbind.data.frame(GPP = sim_list, 
                                        light = sim_pars$X[,2],
                                        discharge = sim_pars$X[,3]))  

  
  # Make the model formula and priors
  bform <- brms::bf(GPP | mi() ~ light + discharge + ar(p = 1))
  bprior <- c(prior(normal(0,1), class = 'ar'),
              prior(normal(0,5), class = 'b'))
  
  oopts <- options(future.globals.maxSize = 1.0 * 1e11)  ## 1e11 bytes = 100 GB (1e10 bytes = 10 GB) - CT
  on.exit(options(oopts))
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
      slice((nrow(dat_full)-forecast_days):nrow(dat_full)) %>%  # Minor point (CT): this gives 74 forecast days (obs 292-365 inclusive); for exactly 73 forecast days, would need to change to (nrow(dat_full) - forecast_days + 1)
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


# fit models
start <- 989
end <- 1000
if (modelOpt) {
  for (i in start:end) {
    CurSim <- i
    #### set up data for this iteration
    
    gauss_sim_CurSim_vals <- sim_list_miss_short[[CurSim]][[3]]
    gauss_sim_CurSim_params <- sim_list_miss_short[[CurSim]][[2]]
    
    # "Full" dataset (i.e. the first element in the gauss_sim_CurSim_vals list)
    # dat_full <- data.frame("date" = 1:365,
    #                        "GPP" = gauss_sim_CurSim_vals$y_noMiss,
    #                        "light" = gauss_sim_CurSim_params$X[,2],
    #                        "discharge" = gauss_sim_CurSim_params$X[,3]
    # )
    
    #####################################################
    #### MODEL RUN KALMAN FILTER ##############
    #########################################################
    arima_kalman_MAR<- fit_arima_Kalman(sim_list = gauss_sim_CurSim_vals,
                                        sim_pars = gauss_sim_CurSim_params, 
                                        forecast = FALSE)
    
    
    ## pull out and label what we need ###
    arimaKalman_MAR_df <- arima_kalman_MAR$arima_pars
    arimaKalman_MAR_df$missingness <- 'MAR'
    arimaKalman_MAR_df$type <- 'Kalman Filter'
    arimaKalman_MAR_df$curSim <- CurSim
    
    #####################################################
    #### MODEL RUN Simple Drop Missing ##############
    #########################################################
    arima_drop_MAR<- fit_arima_dropmissing(sim_list = gauss_sim_CurSim_vals,
                                           sim_pars = gauss_sim_CurSim_params, 
                                           forecast = FALSE)
    
    
    ## pull out and label what we need ###
    arimaDrop_MAR_df <- arima_drop_MAR$arima_pars
    arimaDrop_MAR_df$missingness <- 'MAR'
    arimaDrop_MAR_df$type <- 'simple drop NA'
    arimaDrop_MAR_df$curSim <- CurSim
    
    #####################################################
    #### MODEL RUN Complete Case Drop Missing ##############
    #########################################################
    arima_dropCC_MAR<- fit_arima_dropmissing_CC(sim_list = gauss_sim_CurSim_vals,
                                                sim_pars = gauss_sim_CurSim_params, 
                                                forecast = FALSE)
    
    
    ## pull out and label what we need ###
    arimadropCC_MAR_df <- arima_dropCC_MAR$arima_pars
    arimadropCC_MAR_df$missingness <- 'MAR'
    arimadropCC_MAR_df$type <- 'comple case drop NA'
    arimadropCC_MAR_df$curSim <- CurSim
    
    #####################################################
    #### MODEL RUN Multiple Imputation ##############
    #########################################################
    if (i %in% c(1:100)) {
      arimaMI_MAR_df <- arimadropCC_MAR_df
      arimaMI_MAR_df[, c("param_value", "param_se")] <- NA
      arimaMI_MAR_df$type <- 'Multiple Imputation'
    } else {
      arima_MI_MAR<- fit_arima_MI(sim_list = gauss_sim_CurSim_vals,
                                  sim_pars = gauss_sim_CurSim_params, 
                                  forecast = FALSE,
                                  imputationsnum=5)
      
      ## pull out and label what we need ###
      arimaMI_MAR_df <- arima_MI_MAR[[1]]
      arimaMI_MAR_df$missingness <- 'MAR'
      arimaMI_MAR_df$type <- 'Multiple Imputation'
      arimaMI_MAR_df$curSim <- CurSim
    }
    
    #####################################################
    #### Data Augmentation with brms ##############
    #########################################################
    brms_MAR<- fit_brms_model(sim_list = gauss_sim_CurSim_vals,
                                                sim_pars = gauss_sim_CurSim_params, 
                                                forecast = FALSE)
    
    ## pull out and label what we need ###
    brms_MAR_df <- brms_MAR$brms_pars[[1]]
    brms_MAR_df <- brms_MAR_df %>% 
      filter(parameter != "Intercept")
    brms_MAR_df$missingness <- 'MAR'
    brms_MAR_df$type <- 'data augmentation'
    brms_MAR_df$curSim <- CurSim
    brms_MAR_df <- brms_MAR_df %>% 
      mutate(missingprop_autocor = NA) %>% 
      rename(parameters = parameter,
             param_value = mean,
             param_se = sd
             ) %>% 
      select(names(arimaMI_MAR_df))
    
    ###################################################
    #### COMBINE ALL MODEL RUNS AND SAVE #########
    #################################################
    
    ###
    # put all model parameters together
    params_MAR_all <- rbind(
      arimaKalman_MAR_df, 
      arimaDrop_MAR_df,
      arimadropCC_MAR_df,
      arimaMI_MAR_df,
      brms_MAR_df
    )
    
    # add in the "simulation number" for this iteration (which is stored in the name of the data list element)
    simName <- names(sim_list_miss_short)[CurSim]
    
    params_MAR_all$sim_no <- simName
    
    # add simulation parameters
    simParams_i <- data.frame("simParams" = c(gauss_sim_CurSim_params$phi, gauss_sim_CurSim_params$beta[1], 
                                  gauss_sim_CurSim_params$beta[2], 
                                  gauss_sim_CurSim_params$beta[3], NA),
                              "parameters" = c("phi", "intercept", "xreg1", "xreg2", "sigma"))
    params_MAR_all <- params_MAR_all %>% 
      mutate(parameters = replace(parameters, list = parameters %in% c("ar1", "phi"), "phi"),
             parameters = replace(parameters, list = parameters %in% c("xreg1", "b_light"), "xreg1"),
             parameters = replace(parameters, list = parameters %in% c("xreg2", "b_discharge"), "xreg2"),
             parameters = replace(parameters, list = parameters %in% c("intercept", "b_Intercept"), "intercept"))
    params_MAR_all <- params_MAR_all %>% 
      left_join(simParams_i)
    
    # save name of missing values list (contains info about missingness and autocorrelation)
    params_MAR_all$missingnessName <- names(sim_list_miss_short[[CurSim]][3])
    # save output 
    if(i == start#1
       ) {
      mod_df <- params_MAR_all
    } else {
      mod_df <- mod_df %>% rbind(params_MAR_all)
    }
    # save output
    write.csv(mod_df, paste0("./data/model_results/sampleSizeAnalysis/modelOutput_iteration",start,"_through_",i,".csv"))
  }
} 


# read in data ------------------------------------------------------------

fileNames <- list.files("./data/model_results/sampleSizeAnalysis/")

for (i in seq_along(fileNames)) {
 temp_file <- read.csv(paste0("./data/model_results/sampleSizeAnalysis/",fileNames[i])) 
 if (i == 1) {
   mod_df <- temp_file
 } else {
   mod_df <- mod_df %>% rbind(temp_file)
 }
}

# Calculate parameter recovery and format data ---------------------------------
# calculate parameter recovery
mod_df <- mod_df %>% 
 dplyr::mutate(paramDiff = ((param_value - simParams)/abs(simParams)),
               paramDiff_absDiff = abs((param_value - simParams)/abs(simParams)))
 
# retrieve info about % missingness and initial dataset size
mod_df$amountMiss <- str_extract(mod_df$missingnessName, pattern = "propMissAct_[0-9.]+") %>% 
  str_extract(pattern = "[0-9.]+") %>% 
  as.numeric()
 
# retrieve info about dataset length prior to simulation
mod_df$N <- str_extract(mod_df$sim_no, pattern = "N_[0-9]+") %>% 
  str_extract(pattern = "[0-9]+") %>% 
  as.numeric()

# update names of parameters 
mod_df <- mod_df %>% 
  mutate(parameters = replace(parameters, list = parameters %in% c("ar1", "phi"), "Phi"),
         parameters = replace(parameters, list = parameters %in% c("xreg1", "xreg2"), "Beta covariates"))
# save data 
saveRDS(mod_df,
        file = here::here("data/gauss_ar1_0miss_datasets_ModelResults_forSampleSizeSupplement.rds")
)

# calculate median values per bin of missingness --------------------------------------------
mod_df_bin <- mod_df %>% 
  mutate(as.factor(amountMiss)) %>% 
  group_by(parameters, type, amountMiss) %>% 
  summarize(
    paramDiffAbsDiff_mean = mean(paramDiff_absDiff, na.rm = TRUE),
    paramDiffAbsDiff_med = median(paramDiff_absDiff, na.rm = TRUE),
    paramDiffAbsDiff_SD = sd(paramDiff_absDiff, na.rm = TRUE),
    n_paramDiffAbsDiff = length(paramDiff_absDiff),
    paramDiff_mean = mean(paramDiff, na.rm = TRUE),
    paramDiff_med = median(paramDiff, na.rm = TRUE),
    paramDiff_SD = sd(paramDiff, na.rm = TRUE),
    n_paramDiff = length(paramDiff)#,
    #SE_mean = mean(param_se, na.rm = TRUE)
  ) %>% 
  filter(!parameters %in% c("intercept", "sigma"))

# Figure for median error ------------------------------------------------------------
(gauss_paramRecovery_bias <- ggplot(data = mod_df_bin %>% filter(parameters != "sigma"), aes(x = #as.factor
                                                                                           (amountMiss), y = paramDiff_med)) +
   ggh4x::facet_grid2(rows = vars(factor(parameters, levels = c("Phi", "Beta covariates"))), scales = "free_y", )+
   geom_hline(aes(yintercept = 0), colour = "grey") + 
   #geom_errorbar(aes(ymin=paramDiff_mean - paramDiff_SD, ymax=paramDiff_mean + paramDiff_SD, color = as.factor(type)), 
   #size=0.3, width=0, position = position_dodge(width=0.03))+
   #geom_ribbon(aes(ymin = paramDiff_mean - paramDiff_SD, ymax = paramDiff_mean + paramDiff_SD, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
   #geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
   geom_line(aes(col = type)) +
   geom_point(aes(col = type), alpha = .8#, position = position_dodge(width=0.05)
              ) +
   #geom_boxplot(aes(col = type)) +
   theme_classic() +
   xlab("Proportion of missing data")+ 
   theme(legend.position="top")+
   #theme(legend.title=element_blank())+
   ylab("Median Error")+ 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
  #ylim(c(-1, 1)) #+
   #xlim(c(-0.03,0.65))# + 
   #ylim(c(0,.3)) +
   scale_color_discrete(type = c( "#D55E00",  "#CC79A7", "#0072B2",
                                  "#009E73", "#E69F00"),
    labels = c( "Data Deletion-Complete","Data Augmentation", "Kalman Filter","Multiple Imputation",  
                "Data Deletion-Simple")) + 
   guides(color = guide_legend(title = NULL))
    )


# Figure for Absolute Median Error ----------------------------------------
(gauss_paramRecovery_SE <- ggplot(data = mod_df_bin %>% filter(parameters != "intercept"), aes(x = amountMiss, y = paramDiffAbsDiff_med )) +
   ggh4x::facet_grid2(rows = vars(factor(parameters, levels = c("Phi", "Beta covariates"))), scales = "free_y", )+
   geom_hline(aes(yintercept = 0), colour = "grey") + 
    #geom_errorbar(aes(ymin=paramDiff_mean - paramDiff_SD, ymax=paramDiff_mean + paramDiff_SD, color = as.factor(type)), 
    #size=0.3, width=0, position = position_dodge(width=0.03))+
    #geom_ribbon(aes(ymin = paramDiff_mean - paramDiff_SD, ymax = paramDiff_mean + paramDiff_SD, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
    geom_line(aes(col = type)) + 
    geom_point(aes(col = type)) +
    #geom_boxplot(aes(col = type)) +
    theme_classic() +
    xlab("Proportion of missing data") + 
    theme(legend.position="top")+
    #theme(legend.title=element_blank())+
    ylab("Absolute Median Error")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))+
   scale_color_discrete(type = c( "#D55E00",  "#CC79A7", "#0072B2",
                                  "#009E73", "#E69F00"),
                        labels = c( "Data Deletion-Complete","Data Augmentation", "Kalman Filter","Multiple Imputation",  
                                    "Data Deletion-Simple")) + 
   guides(color = guide_legend(title = NULL))
    #ylim(c(0, 10)) #+
  #xlim(c(-0.03,0.65))# + 
  #ylim(c(0,.3)) +
  #scale_color_discrete(type = c( "#CC79A7",  "#D55E00","#E69F00","#0072B2", "#009E73"),
  #  labels = c( "Data Augmentation",  "Data Deletion-Complete","Data Deletion-Simple","Kalman Filter","Multiple Imputation"))
  # scale_color_discrete(type = c("#E69F00", "#D55E00","#009E73","#0072B2", "#CC79A7"),
  #                      labels = c("Data Deletion-Simple", "Data Deletion-Complete","Multiple Imputation","Kalman Filter", "Data Augmentation"))
)
 

# Figure for Coverage -----------------------------------------------------
mod_df_cov_temp <- mod_df %>% 
  mutate(CI95_lower = param_value - 1.96*param_se,
         CI95_upper = param_value + 1.96*param_se) %>% 
  filter(parameters != "sigma") %>% 
  filter(!is.na(param_se))# randomly there are some model runs that don't have SE? 

# is the true parameter within the 95% CI? 
mod_df_cov_temp$coverage <- c(mod_df_cov_temp$simParams >= mod_df_cov_temp$CI95_lower & 
                                mod_df_cov_temp$simParams <= mod_df_cov_temp$CI95_upper)

## count the # of models w/ and without coverage for each bin of missingness and autocorrelation
mod_df_cov <- mod_df_cov_temp %>% 
  group_by(amountMiss, type, parameters) %>% 
  summarize(coverageNumber = sum(coverage), # the number of models that have coverage
            modelRunN = length(!is.na(coverage))# the total number of models 
  ) %>% 
  mutate(coveragePerc = coverageNumber/modelRunN)


(gauss_paramRecovery_coverage <- ggplot(data = mod_df_cov %>% filter(parameters != "intercept"), aes(x = amountMiss, y = coveragePerc)) +
    #facet_grid(rows = vars(parameters), scale = "free") +
    #geom_col(aes(x = amtMiss, y = coveragePerc, color = as.factor(type), fill = as.factor(type)), 
    # position = "dodge", alpha = .5) + 
    ggh4x::facet_grid2(rows = vars(factor(parameters, levels = c("Phi", "Beta covariates"))), scales = "free_y", )+
    geom_hline(aes(yintercept = .95), colour = "grey", linetype="dotted") +
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("Coverage")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))+
    scale_color_discrete(type = c( "#D55E00",  "#CC79A7", "#0072B2",
                                   "#009E73", "#E69F00"),
                         labels = c( "Data Deletion-Complete","Data Augmentation", "Kalman Filter","Multiple Imputation",  
                                     "Data Deletion-Simple")) +  
    guides(color = guide_legend(title = NULL))
   #xlim(c(-0.03,0.65)) + 
    #ylim(c(0,.3)) +
    # scale_color_discrete(type = c( "#CC79A7",  "#D55E00","#E69F00","#0072B2", "#009E73"),
    #                      labels = c( "Data Augmentation",  "Data Deletion-Complete","Data Deletion-Simple","Kalman Filter","Multiple Imputation"))
    # 
)


# put plots together and save  --------------------------------------------

library(ggpubr)
library(patchwork)
# remove axis labels for first two plots

gauss_paramRecovery_bias2<-gauss_paramRecovery_bias+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), legend.text = element_text(size=7))

gauss_paramRecovery_SE2<-gauss_paramRecovery_SE+ 
  theme(strip.text.x = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),legend.text = element_text(size=7)) +
    guides(color = "none")

gauss_paramRecovery_coverage2<-gauss_paramRecovery_coverage+ theme(strip.text.x = element_blank(),legend.text = element_text(size=7)) +
  guides(color = "none")

sampleSizeAnalysisFig <- gauss_paramRecovery_bias2/gauss_paramRecovery_SE2/gauss_paramRecovery_coverage2 + plot_annotation(tag_levels = c('A'))

png(file = "./figures/sampleSizeAnalysisFigures.png", width = 6.5, height = 8, units = "in", res = 700)
sampleSizeAnalysisFig
dev.off()
