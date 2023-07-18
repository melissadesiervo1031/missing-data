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
# MISSING AT RANDOM low autocor
############################################################################################################


gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_10.rds")
GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]

arima_kalman_MAR<- fit_arima_Kalman(GPP_sim_MAR$y,GPP_sim_MAR$sim_params)


############ formatting for figure #############

names(arima_kalman_MAR) <- names(GPP_sim_MAR[["y"]])

modelNAparamlist<-purrr::map(arima_kalman_MAR , ~.["arima_pars"])
modelNASElist<-purrr::map(arima_kalman_MAR , ~.["arima_errors"])

modelNAparamlist2 <- lapply(modelNAparamlist, function(x) as.data.frame(do.call(rbind, x)))
modelNASElist2 <- lapply(modelNASElist, function(x) as.data.frame(do.call(rbind, x)))


modelNAparamdf <- map_df(modelNAparamlist2, ~as.data.frame(.x), .id="missingprop_autocor")
modelNASEdf <- map_df(modelNASElist2, ~as.data.frame(.x), .id="missingprop_autocor")


modelNAdf<-modelNAparamdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge) %>% mutate(missingness="MAR low auto") %>% mutate(type="Kalman filter")

modelNASEdf<-modelNASEdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge) %>% mutate(missingness="MAR low auto") %>% mutate(type="Kalman filter")

##pull out missing prop and autocorr values ##

modelNAdf$missingprop<-as.numeric(str_extract(modelNAdf$missingprop1, "\\d+\\.*\\d*"))

modelNAdf$autocorr<-as.numeric(str_extract(modelNAdf$autoCorr1, "\\d+\\.*\\d*"))

modelNASEdf$missingprop<-as.numeric(str_extract(modelNASEdf$missingprop1, "\\d+\\.*\\d*"))

modelNASEdf$autocorr<-as.numeric(str_extract(modelNASEdf$autoCorr1, "\\d+\\.*\\d*"))


## long form ##

paramNAlong <- gather(modelNAdf, param, value, ar1:discharge, factor_key=TRUE)

paramNASElong <- gather(modelNASEdf, param, SE, ar1:discharge, factor_key=TRUE)

paramNAlong2<-merge(paramNAlong, paramNASElong)



################################################## #########################################################
# MISSING AT RANDOM low autocor
############################################################################################################


gauss_sim_MAR_datasets_high <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_90.rds")
GPP_sim_MAR_high<- gauss_sim_MAR_datasets_high [[1]]

arima_kalman_MAR_high<- fit_arima_Kalman(GPP_sim_MAR_high$y,GPP_sim_MAR_high$sim_params)


############ formatting for figure #############

names(arima_kalman_MAR_high) <- names(GPP_sim_MAR_high[["y"]])

modelNAparamlist_high<-purrr::map(arima_kalman_MAR_high , ~.["arima_pars"])
modelNASElist_high<-purrr::map(arima_kalman_MAR_high , ~.["arima_errors"])

modelNAparamlist_high2 <- lapply(modelNAparamlist_high, function(x) as.data.frame(do.call(rbind, x)))
modelNASElist_high2 <- lapply(modelNASElist_high, function(x) as.data.frame(do.call(rbind, x)))


modelNAparamdf_high <- map_df(modelNAparamlist_high2, ~as.data.frame(.x), .id="missingprop_autocor")
modelNASEdf_high <- map_df(modelNASElist_high2, ~as.data.frame(.x), .id="missingprop_autocor")


modelNAdf_high<-modelNAparamdf_high  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge) %>% mutate(missingness="MAR high auto") %>% mutate(type="Kalman filter")

modelNASEdf_high<-modelNASEdf_high  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge)%>% mutate(missingness="MAR high auto")  %>% mutate(type="Kalman filter")

##pull out missing prop and autocorr values ##

modelNAdf_high$missingprop<-as.numeric(str_extract(modelNAdf_high$missingprop1, "\\d+\\.*\\d*"))

modelNAdf_high$autocorr<-as.numeric(str_extract(modelNAdf_high$autoCorr1, "\\d+\\.*\\d*"))

modelNASEdf_high$missingprop<-as.numeric(str_extract(modelNASEdf_high$missingprop1, "\\d+\\.*\\d*"))

modelNASEdf_high$autocorr<-as.numeric(str_extract(modelNASEdf_high$autoCorr1, "\\d+\\.*\\d*"))


## long form ##

paramNAlong_high <- gather(modelNAdf_high, param, value, ar1:discharge, factor_key=TRUE)

paramNASElong_high <- gather(modelNASEdf_high, param, SE, ar1:discharge, factor_key=TRUE)

paramNAlong_high2<-merge(paramNAlong_high, paramNASElong_high)




############ MISSING NOT AT RANDOM MNAR ##############################

gauss_sim_MNAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_minMaxMiss.rds")
GPP_sim_MNAR<- gauss_sim_MNAR_datasets [[1]]

arima_kalman_MNAR<- fit_arima_Kalman(GPP_sim_MNAR$y,GPP_sim_MNAR$sim_params)

########### formatting for figure #############

names(arima_kalman_MNAR) <- names(GPP_sim_MAR[["y"]])

modelNAparamlistMNAR<-purrr::map(arima_kalman_MNAR , ~.["arima_pars"])
modelNASElistMNAR<-purrr::map(arima_kalman_MNAR , ~.["arima_errors"])

modelNAparamlistMNAR2 <- lapply(modelNAparamlistMNAR, function(x) as.data.frame(do.call(rbind, x)))
modelNASElistMNAR2 <- lapply(modelNASElistMNAR, function(x) as.data.frame(do.call(rbind, x)))


modelNAparamdfMNAR <- map_df(modelNAparamlistMNAR2, ~as.data.frame(.x), .id="missingprop_autocor")
modelNASEdfMNAR <- map_df(modelNASElistMNAR2, ~as.data.frame(.x), .id="missingprop_autocor")

modelNAdfMNAR<-modelNAparamdfMNAR  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2)  %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge) %>% mutate(missingness="MNAR") %>% mutate(type="Kalman filter")

modelNASEdfMNAR<-modelNASEdfMNAR  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge)%>% mutate(missingness="MNAR") %>% mutate(type="Kalman filter")


##pull out missing prop and autocorr values ##

modelNAdfMNAR$missingprop<-as.numeric(str_extract(modelNAdfMNAR$missingprop1, "\\d+\\.*\\d*"))

modelNAdfMNAR$autocorr<-as.numeric(str_extract(modelNAdfMNAR$autoCorr1, "\\d+\\.*\\d*"))

modelNASEdfMNAR$missingprop<-as.numeric(str_extract(modelNASEdfMNAR$missingprop1, "\\d+\\.*\\d*"))

modelNASEdfMNAR$autocorr<-as.numeric(str_extract(modelNASEdfMNAR$autoCorr1, "\\d+\\.*\\d*"))


# long form #

paramNAlongMNAR <- gather(modelNAdfMNAR, param, value, ar1:discharge, factor_key=TRUE)

paramNASElongMNAR <- gather(modelNASEdfMNAR, param, SE, ar1:discharge, factor_key=TRUE)

paramNAlongMNAR2<-merge(paramNAlongMNAR, paramNASElongMNAR)
