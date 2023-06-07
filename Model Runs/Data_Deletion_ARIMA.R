# Load packages
library(here)
library(tidyverse)
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

arima_drop_MAR<- fit_arima_dropmissing(GPP_sim_MAR$y,GPP_sim_MAR$sim_params)

########### formatting for figure #############

names(arima_drop_MAR) <- names(GPP_sim_MAR)

modeldropparamlist<-purrr::map(arima_drop_MAR , ~.["arima_pars"])
modeldropSElist<-purrr::map(arima_drop_MAR , ~.["arima_errors"])

modeldropparamlist2 <- lapply(modeldropparamlist, function(x) as.data.frame(do.call(rbind, x)))
modeldropSElist2 <- lapply(modeldropSElist, function(x) as.data.frame(do.call(rbind, x)))


modeldropparamdf <- map_df(modeldropparamlist2, ~as.data.frame(.x), .id="missingprop")
modeldropSEdf <- map_df(modeldropSElist2, ~as.data.frame(.x), .id="missingprop")

modeldropdf<-modeldropparamdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Data Deletion")

modeldropSEdf<-modeldropSEdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Data Deletion")


missingprop<-seq(0.00, 0.95, by=0.05)

modeldropdf2<-cbind(missingprop, modeldropdf)

SEdropdf2<-cbind(missingprop, modeldropSEdf)


paramdroplong <- gather(modeldropdf2, param, value, ar1:discharge, factor_key=TRUE)

paramdropSElong <- gather(SEdropdf2, param, SE, ar1:discharge, factor_key=TRUE)

paramdroplong2<-merge(paramdroplong, paramdropSElong)

########################################################
############ MISSING NOT AT RANDOM MNAR ##############################
####################################################################



gauss_sim_MNAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_minMaxMiss.rds")
GPP_sim_MNAR<- gauss_sim_MNAR_datasets [[1]]

arima_drop_MNAR<- fit_arima_dropmissing(GPP_sim_MNAR$y,GPP_sim_MNAR$sim_params)

########### formatting for figure #############

names(arima_drop_MNAR) <- names(GPP_sim_MAR)

modeldropparamlistMNAR<-purrr::map(arima_drop_MNAR , ~.["arima_pars"])
modeldropSElistMNAR<-purrr::map(arima_drop_MNAR , ~.["arima_errors"])

modeldropparamlistMNAR2 <- lapply(modeldropparamlistMNAR, function(x) as.data.frame(do.call(rbind, x)))
modeldropSElistMNAR2 <- lapply(modeldropSElistMNAR, function(x) as.data.frame(do.call(rbind, x)))


modeldropparamdfMNAR <- map_df(modeldropparamlistMNAR2, ~as.data.frame(.x), .id="missingprop")
modeldropSEdfMNAR <- map_df(modeldropSElistMNAR2, ~as.data.frame(.x), .id="missingprop")

modeldropdfMNAR<-modeldropparamdfMNAR  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Data Deletion")

modeldropSEdfMNAR<-modeldropSEdfMNAR  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Data Deletion")


missingprop<-seq(0.00, 0.95, by=0.05)

modeldropdfMNAR2<-cbind(missingprop, modeldropdfMNAR)

SEdropdfMNAR2<-cbind(missingprop, modeldropSEdfMNAR)


paramdroplongMNAR <- gather(modeldropdfMNAR2, param, value, ar1:discharge, factor_key=TRUE)

paramdropSElongMNAR <- gather(SEdropdfMNAR2, param, SE, ar1:discharge, factor_key=TRUE)

paramdroplongMNAR2<-merge(paramdroplongMNAR, paramdropSElongMNAR)

