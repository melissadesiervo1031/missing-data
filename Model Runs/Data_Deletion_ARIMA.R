# Load packages
library(here)
library(tidyverse)
library(forecast)
library(xts)
library(nlme)
library(tidyverse)
library(lubridate)



################################################## #########################################################
# MISSING AT RANDOM low autocor ##
############################################################################################################

gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_10.rds")
GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]

arima_drop_MAR<- fit_arima_dropmissing(GPP_sim_MAR$y,GPP_sim_MAR$sim_params)

########### formatting for figure #############

names(arima_drop_MAR) <- names(GPP_sim_MAR[["y"]])

modeldropparamlist<-purrr::map(arima_drop_MAR , ~.["arima_pars"])
modeldropSElist<-purrr::map(arima_drop_MAR , ~.["arima_errors"])

modeldropparamlist2 <- lapply(modeldropparamlist, function(x) as.data.frame(do.call(rbind, x)))
modeldropSElist2 <- lapply(modeldropSElist, function(x) as.data.frame(do.call(rbind, x)))


modeldropparamdf <- map_df(modeldropparamlist2, ~as.data.frame(.x), .id="missingprop_autocor")
modeldropSEdf <- map_df(modeldropSElist2, ~as.data.frame(.x), .id="missingprop_autocor")

modeldropdf<-modeldropparamdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge)  %>% mutate(missingness="MAR low auto") %>% mutate(type="Data Deletion")

modeldropSEdf<-modeldropSEdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2)%>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge) %>% mutate(missingness="MAR low auto") %>% mutate(type="Data Deletion")

## pull out missing prop and autocorr values ##

modeldropdf$missingprop<-as.numeric(str_extract(modeldropdf$missingprop1, "\\d+\\.*\\d*"))

modeldropdf$autocorr<-as.numeric(str_extract(modeldropdf$autoCorr1, "\\d+\\.*\\d*"))

modeldropSEdf$missingprop<-as.numeric(str_extract(modeldropSEdf$missingprop1, "\\d+\\.*\\d*"))

modeldropSEdf$autocorr<-as.numeric(str_extract(modeldropSEdf$autoCorr1, "\\d+\\.*\\d*"))

## long form ##

paramdroplong <- gather(modeldropdf, param, value, ar1:discharge, factor_key=TRUE)

paramdropSElong <- gather(modeldropSEdf, param, SE, ar1:discharge, factor_key=TRUE)

paramdroplong2<-merge(paramdroplong, paramdropSElong)



################################################## #########################################################
# MISSING AT RANDOM HIGH autocor ##
############################################################################################################

gauss_sim_MAR_high_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_90.rds")
GPP_sim_MAR_high<- gauss_sim_MAR_high_datasets [[1]]

arima_drop_MAR_high<- fit_arima_dropmissing(GPP_sim_MAR_high$y,GPP_sim_MAR_high$sim_params)

########### formatting for figure #############

names(arima_drop_MAR_high) <- names(GPP_sim_MAR_high[["y"]])

modeldropparamlist_high<-purrr::map(arima_drop_MAR_high , ~.["arima_pars"])
modeldropSElist_high<-purrr::map(arima_drop_MAR_high , ~.["arima_errors"])

modeldropparamlist_high2 <- lapply(modeldropparamlist_high, function(x) as.data.frame(do.call(rbind, x)))
modeldropSElist_high2 <- lapply(modeldropSElist_high, function(x) as.data.frame(do.call(rbind, x)))


modeldropparamdf_high <- map_df(modeldropparamlist_high2, ~as.data.frame(.x), .id="missingprop_autocor")
modeldropSEdf_high <- map_df(modeldropSElist_high2, ~as.data.frame(.x), .id="missingprop_autocor")

modeldropdf_high<-modeldropparamdf_high  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge) %>% mutate(missingness="MAR high auto")%>% mutate(type="Data Deletion")

modeldropSEdf_high<-modeldropSEdf_high  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2)%>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge) %>% mutate(missingness="MAR high auto") %>% mutate(type="Data Deletion")

## pull out missing prop and autocorr values ##

modeldropdf_high$missingprop<-as.numeric(str_extract(modeldropdf_high$missingprop1, "\\d+\\.*\\d*"))

modeldropdf_high$autocorr<-as.numeric(str_extract(modeldropdf_high$autoCorr1, "\\d+\\.*\\d*"))

modeldropSEdf_high$missingprop<-as.numeric(str_extract(modeldropSEdf_high$missingprop1, "\\d+\\.*\\d*"))

modeldropSEdf_high$autocorr<-as.numeric(str_extract(modeldropSEdf_high$autoCorr1, "\\d+\\.*\\d*"))

## long form ##

paramdroplong_high <- gather(modeldropdf_high, param, value, ar1:discharge, factor_key=TRUE)

paramdropSElong_high <- gather(modeldropSEdf_high, param, SE, ar1:discharge, factor_key=TRUE)

paramdroplong_high2<-merge(paramdroplong_high, paramdropSElong_high)





########################################################
############ MISSING NOT AT RANDOM MNAR ##############################
####################################################################



gauss_sim_MNAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_minMaxMiss.rds")
GPP_sim_MNAR<- gauss_sim_MNAR_datasets [[1]]

arima_drop_MNAR<- fit_arima_dropmissing(GPP_sim_MNAR$y,GPP_sim_MNAR$sim_params)

########### formatting for figure #############

names(arima_drop_MNAR) <- names(GPP_sim_MAR[["y"]])

modeldropparamlistMNAR<-purrr::map(arima_drop_MNAR , ~.["arima_pars"])
modeldropSElistMNAR<-purrr::map(arima_drop_MNAR , ~.["arima_errors"])

modeldropparamlistMNAR2 <- lapply(modeldropparamlistMNAR, function(x) as.data.frame(do.call(rbind, x)))
modeldropSElistMNAR2 <- lapply(modeldropSElistMNAR, function(x) as.data.frame(do.call(rbind, x)))


modeldropparamdfMNAR <- map_df(modeldropparamlistMNAR2, ~as.data.frame(.x), .id="missingprop_autocor")
modeldropSEdfMNAR <- map_df(modeldropSElistMNAR2, ~as.data.frame(.x), .id="missingprop_autocor")

modeldropdfMNAR<-modeldropparamdfMNAR  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2)  %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge)  %>% mutate(missingness="MNAR") %>% mutate(type="Data Deletion")

modeldropSEdfMNAR<-modeldropSEdfMNAR  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>%   separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge)  %>% mutate(missingness="MNAR") %>% mutate(type="Data Deletion")

## pull out missing prop and autocor values ##

modeldropdfMNAR$missingprop<-as.numeric(str_extract(modeldropdfMNAR$missingprop1, "\\d+\\.*\\d*"))

modeldropdfMNAR$autocorr<-as.numeric(str_extract(modeldropdfMNAR$autoCorr1, "\\d+\\.*\\d*"))

modeldropSEdfMNAR$missingprop<-as.numeric(str_extract(modeldropSEdfMNAR$missingprop1, "\\d+\\.*\\d*"))

modeldropSEdfMNAR$autocorr<-as.numeric(str_extract(modeldropSEdfMNAR$autoCorr1, "\\d+\\.*\\d*"))




paramdroplongMNAR <- gather(modeldropdfMNAR, param, value, ar1:discharge, factor_key=TRUE)

paramdropSElongMNAR <- gather(modeldropSEdfMNAR, param, SE, ar1:discharge, factor_key=TRUE)

paramdroplongMNAR2<-merge(paramdroplongMNAR, paramdropSElongMNAR)


### combine all three MAR, MAR, MNAR versions into 1 ###


paramdropall<-rbind(paramdroplong2,paramdroplong_high2, paramdroplongMNAR2)
