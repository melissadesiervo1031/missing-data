# Load packages
library(here)
library(tidyverse)
library(Amelia)
library(forecast)
library(xts)
library(nlme)
library(tidyverse)
library(lubridate)


################################################## #########################################################
# MISSING AT RANDOM
############################################################################################################


gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss.rds")
GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]

arima_kalman_MAR<- fit_arima_Kalman(GPP_sim_MAR$y,GPP_sim_MAR$sim_params)


############ formatting for figure #############

names(arima_kalman_MAR) <- names(GPP_sim_MAR)

modelNAparamlist<-purrr::map(arima_kalman_MAR , ~.["arima_pars"])
modelNASElist<-purrr::map(arima_kalman_MAR , ~.["arima_errors"])

modelNAparamlist2 <- lapply(modelNAparamlist, function(x) as.data.frame(do.call(rbind, x)))
modelNASElist2 <- lapply(modelNASElist, function(x) as.data.frame(do.call(rbind, x)))


modelNAparamdf <- map_df(modelNAparamlist2, ~as.data.frame(.x), .id="missingprop")
modelNASEdf <- map_df(modelNASElist2, ~as.data.frame(.x), .id="missingprop")

modelNAdf<-modelNAparamdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Kalman filter")

modelNASEdf<-modelNASEdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Kalman filter")


missingprop<-seq(0.00, 0.95, by=0.05)

modelNAdf2<-cbind(missingprop, modelNAdf)

SENAdf2<-cbind(missingprop, modelNASEdf)


paramNAlong <- gather(modelNAdf2, param, value, ar1:discharge, factor_key=TRUE)

paramNASElong <- gather(SENAdf2, param, SE, ar1:discharge, factor_key=TRUE)

paramNAlong2<-merge(paramNAlong, paramNASElong)



############ MISSING NOT AT RANDOM MNAR ##############################

gauss_sim_MNAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_minMaxMiss.rds")
GPP_sim_MNAR<- gauss_sim_MNAR_datasets [[1]]

arima_kalman_MNAR<- fit_arima_Kalman(GPP_sim_MNAR$y,GPP_sim_MNAR$sim_params)

########### formatting for figure #############

names(arima_kalman_MNAR) <- names(GPP_sim_MAR)

modelNAparamlistMNAR<-purrr::map(arima_kalman_MNAR , ~.["arima_pars"])
modelNASElistMNAR<-purrr::map(arima_kalman_MNAR , ~.["arima_errors"])

modelNAparamlistMNAR2 <- lapply(modelNAparamlistMNAR, function(x) as.data.frame(do.call(rbind, x)))
modelNASElistMNAR2 <- lapply(modelNASElistMNAR, function(x) as.data.frame(do.call(rbind, x)))


modelNAparamdfMNAR <- map_df(modelNAparamlistMNAR2, ~as.data.frame(.x), .id="missingprop")
modelNASEdfMNAR <- map_df(modelNASElistMNAR2, ~as.data.frame(.x), .id="missingprop")

modelNAdfMNAR<-modelNAparamdfMNAR  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Kalman filter")

modelNASEdfMNAR<-modelNASEdfMNAR  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Kalman filter")


missingprop<-seq(0.00, 0.95, by=0.05)

modelNAdfMNAR2<-cbind(missingprop, modelNAdfMNAR)

SENAdfMNAR2<-cbind(missingprop, modelNASEdfMNAR)


paramNAlongMNAR <- gather(modelNAdfMNAR2, param, value, ar1:discharge, factor_key=TRUE)

paramNASElongMNAR <- gather(SENAdfMNAR2, param, SE, ar1:discharge, factor_key=TRUE)

paramNAlongMNAR2<-merge(paramNAlongMNAR, paramNASElongMNAR)
