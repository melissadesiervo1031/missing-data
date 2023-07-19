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
# MISSING AT RANDOM low auto
############################################################################################################


gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_10.rds")
GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]

arima_mi_MAR<- fit_arima_MI(GPP_sim_MAR$y,GPP_sim_MAR$sim_params, imputationsnum=5)

##pulls out parameters and ses ##

paramlistsim<-arima_mi_MAR[[1]]

selistsim<-arima_mi_MAR[[2]]

avgparamdf <- map_df(paramlistsim, ~as.data.frame(.x), .id="missingprop_autocor")
avglSEdf <- map_df(selistsim, ~as.data.frame(.x), .id="missingprop_autocor")


avgparamdf2<-avgparamdf %>% dplyr::rename(ar1=q.mi.ar1, intercept=q.mi.intercept, light=q.mi.xreg1, discharge=q.mi.xreg2) %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge)  %>% mutate(missingness="MAR low auto") %>% mutate(type="Multiple imputations")

avglSEdf2<-avglSEdf  %>% dplyr::rename(ar1=se.mi.ar1, intercept=se.mi.intercept, light=se.mi.xreg1, discharge=se.mi.xreg2) %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge)   %>% mutate(missingness="MAR low auto") %>% mutate(type="Multiple imputations")


## pull out missing prop and autocorr values ##

avgparamdf2$missingprop<-as.numeric(str_extract(avgparamdf2$missingprop1, "\\d+\\.*\\d*"))

avgparamdf2$autocorr<-as.numeric(str_extract(avgparamdf2$autoCorr1, "\\d+\\.*\\d*"))

avglSEdf2$missingprop<-as.numeric(str_extract(avglSEdf2$missingprop1, "\\d+\\.*\\d*"))

avglSEdf2$autocorr<-as.numeric(str_extract(avglSEdf2$autoCorr1, "\\d+\\.*\\d*"))



paramMIlong <- gather(avgparamdf2, param, value, ar1:discharge, factor_key=TRUE)

paramMISElong <- gather(avglSEdf2, param, SE, ar1:discharge, factor_key=TRUE)

paramMIlong2<-merge(paramMIlong,paramMISElong)


################################################# #########################################################
# MISSING AT RANDOM high auto
############################################################################################################


gauss_sim_MAR_datasets_high <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_90.rds")
GPP_sim_MAR_high<- gauss_sim_MAR_datasets_high [[1]]

arima_mi_MAR_high<- fit_arima_MI(GPP_sim_MAR_high$y,GPP_sim_MAR_high$sim_params, imputationsnum=5)

##pulls out parameters and ses ##

paramlistsim_high<-arima_mi_MAR_high[[1]]

selistsim_high<-arima_mi_MAR_high[[2]]

avgparamdf_high <- map_df(paramlistsim_high, ~as.data.frame(.x), .id="missingprop_autocor")
avglSEdf_high <- map_df(selistsim_high, ~as.data.frame(.x), .id="missingprop_autocor")


avgparamdf2_high<-avgparamdf_high %>% dplyr::rename(ar1=q.mi.ar1, intercept=q.mi.intercept, light=q.mi.xreg1, discharge=q.mi.xreg2) %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge)  %>% mutate(missingness="MAR high auto") %>% mutate(type="Multiple imputations")

avglSEdf2_high<-avglSEdf_high  %>% dplyr::rename(ar1=se.mi.ar1, intercept=se.mi.intercept, light=se.mi.xreg1, discharge=se.mi.xreg2) %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge)  %>% mutate(missingness="MAR high auto") %>% mutate(type="Multiple imputations")


## pull out missing prop and autocorr values ##

avgparamdf2_high$missingprop<-as.numeric(str_extract(avgparamdf2_high$missingprop1, "\\d+\\.*\\d*"))

avgparamdf2_high$autocorr<-as.numeric(str_extract(avgparamdf2_high$autoCorr1, "\\d+\\.*\\d*"))

avglSEdf2_high$missingprop<-as.numeric(str_extract(avglSEdf2_high$missingprop1, "\\d+\\.*\\d*"))

avglSEdf2_high$autocorr<-as.numeric(str_extract(avglSEdf2_high$autoCorr1, "\\d+\\.*\\d*"))



paramMIlong_high <- gather(avgparamdf2_high, param, value, ar1:discharge, factor_key=TRUE)

paramMISElong_high <- gather(avglSEdf2_high, param, SE, ar1:discharge, factor_key=TRUE)

paramMIlong_high2<-merge(paramMIlong_high,paramMISElong_high)


############ MISSING NOT AT RANDOM MNAR ##############################

gauss_sim_MNAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_minMaxMiss.rds")
GPP_sim_MNAR<- gauss_sim_MNAR_datasets [[1]]

arima_mi_MNAR<- fit_arima_MI(GPP_sim_MNAR$y,GPP_sim_MNAR$sim_params, imputationsnum=5)

##pulls out parameters and ses ##

paramlistsimMNAR<-arima_mi_MNAR[[1]]

selistsimMNAR<-arima_mi_MNAR[[2]]


avgparamdfMNAR <- map_df(paramlistsimMNAR, ~as.data.frame(.x), .id="missingprop_autocor")
avglSEdfMNAR <- map_df(selistsimMNAR, ~as.data.frame(.x), .id="missingprop_autocor")


avgparamdfMNAR2<-avgparamdfMNAR %>% dplyr::rename(ar1=q.mi.ar1, intercept=q.mi.intercept, light=q.mi.xreg1, discharge=q.mi.xreg2) %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge)  %>% mutate(missingness="MNAR") %>% mutate(type="Multiple imputations")

avglSEdfMNAR2<-avglSEdfMNAR  %>% dplyr::rename(ar1=se.mi.ar1, intercept=se.mi.intercept, light=se.mi.xreg1, discharge=se.mi.xreg2)%>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) %>% select(missingprop1, autoCorr1, ar1, intercept, light, discharge)  %>% mutate(missingness="MNAR")%>% mutate(type="Multiple imputations")


avgparamdfMNAR2$missingprop<-as.numeric(str_extract(avgparamdfMNAR2$missingprop1, "\\d+\\.*\\d*"))

avgparamdfMNAR2$autocorr<-as.numeric(str_extract(avgparamdfMNAR2$autoCorr1, "\\d+\\.*\\d*"))

avglSEdfMNAR2$missingprop<-as.numeric(str_extract(avglSEdfMNAR2$missingprop1, "\\d+\\.*\\d*"))

avglSEdfMNAR2$autocorr<-as.numeric(str_extract(avglSEdfMNAR2$autoCorr1, "\\d+\\.*\\d*"))



# long form #

paramMIlongMNAR <- gather(avgparamdfMNAR2, param, value, ar1:discharge, factor_key=TRUE)

paramMISElongMNAR <- gather(avglSEdfMNAR2, param, SE, ar1:discharge, factor_key=TRUE)

paramMIlongMNAR2<-merge(paramMIlongMNAR,paramMISElongMNAR)


#combine MAR, MAR, MNAR into one ###

paramMIall<-rbind(paramMIlongMNAR2, paramMIlong_high2, paramMIlong2)