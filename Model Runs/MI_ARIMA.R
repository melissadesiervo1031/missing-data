# Load packages
library(here)
library(tidyverse)
library(Amelia)
library(forecast)
library(xts)
library(nlme)
library(tidyverse)
library(lubridate)

################################################# #########################################################
# MISSING AT RANDOM
############################################################################################################


gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss.rds")
GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]

arima_mi_MAR<- fit_arima_MI(GPP_sim_MAR$y,GPP_sim_MAR$sim_params, imputationsnum=5)

##pulls out parameters and ses ##

paramlistsim<-arima_mi_MAR[[1]]

selistsim<-arima_mi_MAR[[2]]

avgparamdf <- map_df(paramlistsim, ~as.data.frame(.x), .id="missingprop")
avglSEdf <- map_df(selistsim, ~as.data.frame(.x), .id="missingprop")

missingprop<-seq(from=0.00, to =0.95, by=0.05)

avgparamdf2<-avgparamdf %>% dplyr::rename(ar1=q.mi.ar1, intercept=q.mi.intercept, light=q.mi.xreg1, discharge=q.mi.xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Multiple imputations")

avglSEdf2<-avglSEdf  %>% dplyr::rename(ar1=se.mi.ar1, intercept=se.mi.intercept, light=se.mi.xreg1, discharge=se.mi.xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Multiple imputations")


missingprop<-seq(0.00, 0.95, by=0.05)

MIdf2<-cbind(missingprop, avgparamdf2)

MISEdf2<-cbind(missingprop, avglSEdf2)


paramMIlong <- gather(MIdf2, param, value, ar1:discharge, factor_key=TRUE)

paramMISElong <- gather(MISEdf2, param, SE, ar1:discharge, factor_key=TRUE)

paramMIlong2<-merge(paramMIlong,paramMISElong)


############ MISSING NOT AT RANDOM MNAR ##############################

gauss_sim_MNAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_minMaxMiss.rds")
GPP_sim_MNAR<- gauss_sim_MNAR_datasets [[1]]

arima_mi_MNAR<- fit_arima_MI(GPP_sim_MNAR$y,GPP_sim_MNAR$sim_params, imputationsnum=5)

##pulls out parameters and ses ##

paramlistsimMNAR<-arima_mi_MNAR[[1]]

selistsimMNAR<-arima_mi_MNAR[[2]]


avgparamdfMNAR <- map_df(paramlistsimMNAR, ~as.data.frame(.x), .id="missingprop")
avglSEdfMNAR <- map_df(selistsimMNAR, ~as.data.frame(.x), .id="missingprop")

missingprop<-seq(from=0.00, to =0.95, by=0.05)

avgparamdfMNAR2<-avgparamdfMNAR %>% dplyr::rename(ar1=q.mi.ar1, intercept=q.mi.intercept, light=q.mi.xreg1, discharge=q.mi.xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Multiple imputations")

avglSEdfMNAR2<-avglSEdfMNAR  %>% dplyr::rename(ar1=se.mi.ar1, intercept=se.mi.intercept, light=se.mi.xreg1, discharge=se.mi.xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Multiple imputations")


missingprop<-seq(0.00, 0.95, by=0.05)

MIdfMNAR2<-cbind(missingprop, avgparamdfMNAR2)

MISEdfMNAR2<-cbind(missingprop, avglSEdfMNAR2)


paramMIlongMNAR <- gather(MIdfMNAR2, param, value, ar1:discharge, factor_key=TRUE)

paramMISElongMNAR <- gather(MISEdfMNAR2, param, SE, ar1:discharge, factor_key=TRUE)

paramMIlongMNAR2<-merge(paramMIlongMNAR,paramMISElongMNAR)

