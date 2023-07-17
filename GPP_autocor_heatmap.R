## draft of the GPP heat map figure ##

# Load packages
library(here)
library(tidyverse)
library(rstan)
library(Amelia)
library(forecast)
library(xts)
library(nlme)
library(tidyverse)
library(lubridate)
library(here)
library(viridis)

#####Source code for the ARIMA functions, just Kalman filter for now ###

source(here::here("Functions/Arima_Kalman_function.R"))

### Read in MAR datasets with different levels of autocorrellation ##

gauss_sim_MAR_cor_01_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_01.rds")
gauss_sim_MAR_cor_10_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_10.rds")
gauss_sim_MAR_cor_20_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_20.rds")
gauss_sim_MAR_cor_30_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_30.rds")
gauss_sim_MAR_cor_40_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_40.rds")
gauss_sim_MAR_cor_50_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_50.rds")
gauss_sim_MAR_cor_60_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_60.rds")
gauss_sim_MAR_cor_70_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_70.rds")
gauss_sim_MAR_cor_80_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_80.rds")
gauss_sim_MAR_cor_90_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_90.rds")


##pull out the first of each one ##

gauss_sim_MAR_cor_01_1<- gauss_sim_MAR_cor_01_datasets [[1]]
gauss_sim_MAR_cor_10_1<- gauss_sim_MAR_cor_10_datasets [[1]]
gauss_sim_MAR_cor_20_1<- gauss_sim_MAR_cor_20_datasets [[1]]
gauss_sim_MAR_cor_30_1<- gauss_sim_MAR_cor_30_datasets [[1]]
gauss_sim_MAR_cor_40_1<- gauss_sim_MAR_cor_40_datasets [[1]]
gauss_sim_MAR_cor_50_1<- gauss_sim_MAR_cor_50_datasets [[1]]
gauss_sim_MAR_cor_60_1<- gauss_sim_MAR_cor_60_datasets [[1]]
gauss_sim_MAR_cor_70_1<- gauss_sim_MAR_cor_70_datasets [[1]]
gauss_sim_MAR_cor_80_1<- gauss_sim_MAR_cor_80_datasets [[1]]
gauss_sim_MAR_cor_90_1<- gauss_sim_MAR_cor_90_datasets [[1]]

### true estimates without missing data ##

## no missing###


realphi<-gauss_sim_MAR_cor_01_1[["sim_params"]][["phi"]]
realbetas<-gauss_sim_MAR_cor_01_1[["sim_params"]][["beta"]]

trueestdf <- data.frame (param = c("ar1", "intercept", "light", "discharge"), realvalue = c(realphi, realbetas))


## Run the Kalman filter function over each of the autocorrellated missing datasets##

arima_kalman_MAR_cor_01<- fit_arima_Kalman(gauss_sim_MAR_cor_01_1$y,gauss_sim_MAR_cor_01_1$sim_params)
names(arima_kalman_MAR_cor_01) <- names(gauss_sim_MAR_cor_01_1[["y"]]) #names for bookeeping#

arima_kalman_MAR_cor_10<- fit_arima_Kalman(gauss_sim_MAR_cor_10_1$y,gauss_sim_MAR_cor_10_1$sim_params)
names(arima_kalman_MAR_cor_10) <- names(gauss_sim_MAR_cor_10_1[["y"]]) #names for bookeeping#

arima_kalman_MAR_cor_20<- fit_arima_Kalman(gauss_sim_MAR_cor_20_1$y,gauss_sim_MAR_cor_20_1$sim_params)
names(arima_kalman_MAR_cor_20) <- names(gauss_sim_MAR_cor_20_1[["y"]]) #names for bookeeping#

arima_kalman_MAR_cor_30<- fit_arima_Kalman(gauss_sim_MAR_cor_30_1$y,gauss_sim_MAR_cor_30_1$sim_params)
names(arima_kalman_MAR_cor_30) <- names(gauss_sim_MAR_cor_30_1[["y"]]) #names for bookeeping#

arima_kalman_MAR_cor_40<- fit_arima_Kalman(gauss_sim_MAR_cor_40_1$y,gauss_sim_MAR_cor_40_1$sim_params)
names(arima_kalman_MAR_cor_40) <- names(gauss_sim_MAR_cor_40_1[["y"]]) #names for bookeeping#

arima_kalman_MAR_cor_50<- fit_arima_Kalman(gauss_sim_MAR_cor_50_1$y,gauss_sim_MAR_cor_50_1$sim_params)
names(arima_kalman_MAR_cor_50) <- names(gauss_sim_MAR_cor_50_1[["y"]]) #names for bookeeping#


arima_kalman_MAR_cor_60<- fit_arima_Kalman(gauss_sim_MAR_cor_60_1$y,gauss_sim_MAR_cor_60_1$sim_params)
names(arima_kalman_MAR_cor_60) <- names(gauss_sim_MAR_cor_60_1[["y"]]) #names for bookeeping#


arima_kalman_MAR_cor_70<- fit_arima_Kalman(gauss_sim_MAR_cor_70_1$y,gauss_sim_MAR_cor_70_1$sim_params)
names(arima_kalman_MAR_cor_70) <- names(gauss_sim_MAR_cor_70_1[["y"]]) #names for bookeeping#


arima_kalman_MAR_cor_80<- fit_arima_Kalman(gauss_sim_MAR_cor_80_1$y,gauss_sim_MAR_cor_80_1$sim_params)
names(arima_kalman_MAR_cor_80) <- names(gauss_sim_MAR_cor_80_1[["y"]]) #names for bookeeping#


arima_kalman_MAR_cor_90<- fit_arima_Kalman(gauss_sim_MAR_cor_90_1$y,gauss_sim_MAR_cor_90_1$sim_params)
names(arima_kalman_MAR_cor_90) <- names(gauss_sim_MAR_cor_90_1[["y"]]) #names for bookeeping#

## pull out what we need ###

#01#
arima_kalman_MAR_cor_01_param<-purrr::map(arima_kalman_MAR_cor_01 , ~.["arima_pars"])
arima_kalman_MAR_cor_01_error<-purrr::map(arima_kalman_MAR_cor_01 , ~.["arima_errors"])

arima_kalman_MAR_cor_01_paramdf <- lapply(arima_kalman_MAR_cor_01_param, function(x) as.data.frame(do.call(rbind, x)))
arima_kalman_MAR_cor_01_errordf<- lapply(arima_kalman_MAR_cor_01_error, function(x) as.data.frame(do.call(rbind, x)))

kalman_MAR_cor_01_paramdf <- map_df(arima_kalman_MAR_cor_01_paramdf, ~as.data.frame(.x), .id="missingprop_autocor")
kalman_MAR_cor_01_errordf <- map_df(arima_kalman_MAR_cor_01_errordf, ~as.data.frame(.x), .id="missingprop_autocor")

#10#
arima_kalman_MAR_cor_10_param<-purrr::map(arima_kalman_MAR_cor_10 , ~.["arima_pars"])
arima_kalman_MAR_cor_10_error<-purrr::map(arima_kalman_MAR_cor_10 , ~.["arima_errors"])

arima_kalman_MAR_cor_10_paramdf <- lapply(arima_kalman_MAR_cor_10_param, function(x) as.data.frame(do.call(rbind, x)))
arima_kalman_MAR_cor_10_errordf<- lapply(arima_kalman_MAR_cor_10_error, function(x) as.data.frame(do.call(rbind, x)))

kalman_MAR_cor_10_paramdf <- map_df(arima_kalman_MAR_cor_10_paramdf, ~as.data.frame(.x), .id="missingprop_autocor")
kalman_MAR_cor_10_errordf <- map_df(arima_kalman_MAR_cor_10_errordf, ~as.data.frame(.x), .id="missingprop_autocor")


#20#
arima_kalman_MAR_cor_20_param<-purrr::map(arima_kalman_MAR_cor_20 , ~.["arima_pars"])
arima_kalman_MAR_cor_20_error<-purrr::map(arima_kalman_MAR_cor_20 , ~.["arima_errors"])

arima_kalman_MAR_cor_20_paramdf <- lapply(arima_kalman_MAR_cor_20_param, function(x) as.data.frame(do.call(rbind, x)))
arima_kalman_MAR_cor_20_errordf<- lapply(arima_kalman_MAR_cor_20_error, function(x) as.data.frame(do.call(rbind, x)))

kalman_MAR_cor_20_paramdf <- map_df(arima_kalman_MAR_cor_20_paramdf, ~as.data.frame(.x), .id="missingprop_autocor")
kalman_MAR_cor_20_errordf <- map_df(arima_kalman_MAR_cor_20_errordf, ~as.data.frame(.x), .id="missingprop_autocor")


#30#
arima_kalman_MAR_cor_30_param<-purrr::map(arima_kalman_MAR_cor_30 , ~.["arima_pars"])
arima_kalman_MAR_cor_30_error<-purrr::map(arima_kalman_MAR_cor_30 , ~.["arima_errors"])

arima_kalman_MAR_cor_30_paramdf <- lapply(arima_kalman_MAR_cor_30_param, function(x) as.data.frame(do.call(rbind, x)))
arima_kalman_MAR_cor_30_errordf<- lapply(arima_kalman_MAR_cor_30_error, function(x) as.data.frame(do.call(rbind, x)))

kalman_MAR_cor_30_paramdf <- map_df(arima_kalman_MAR_cor_30_paramdf, ~as.data.frame(.x), .id="missingprop_autocor")
kalman_MAR_cor_30_errordf <- map_df(arima_kalman_MAR_cor_30_errordf, ~as.data.frame(.x), .id="missingprop_autocor")


#40#
arima_kalman_MAR_cor_40_param<-purrr::map(arima_kalman_MAR_cor_40 , ~.["arima_pars"])
arima_kalman_MAR_cor_40_error<-purrr::map(arima_kalman_MAR_cor_40 , ~.["arima_errors"])

arima_kalman_MAR_cor_40_paramdf <- lapply(arima_kalman_MAR_cor_40_param, function(x) as.data.frame(do.call(rbind, x)))
arima_kalman_MAR_cor_40_errordf<- lapply(arima_kalman_MAR_cor_40_error, function(x) as.data.frame(do.call(rbind, x)))

kalman_MAR_cor_40_paramdf <- map_df(arima_kalman_MAR_cor_40_paramdf, ~as.data.frame(.x), .id="missingprop_autocor")
kalman_MAR_cor_40_errordf <- map_df(arima_kalman_MAR_cor_40_errordf, ~as.data.frame(.x), .id="missingprop_autocor")


#50$
arima_kalman_MAR_cor_50_param<-purrr::map(arima_kalman_MAR_cor_50 , ~.["arima_pars"])
arima_kalman_MAR_cor_50_error<-purrr::map(arima_kalman_MAR_cor_50 , ~.["arima_errors"])

arima_kalman_MAR_cor_50_paramdf <- lapply(arima_kalman_MAR_cor_50_param, function(x) as.data.frame(do.call(rbind, x)))
arima_kalman_MAR_cor_50_errordf<- lapply(arima_kalman_MAR_cor_50_error, function(x) as.data.frame(do.call(rbind, x)))

kalman_MAR_cor_50_paramdf <- map_df(arima_kalman_MAR_cor_50_paramdf, ~as.data.frame(.x), .id="missingprop_autocor")
kalman_MAR_cor_50_errordf <- map_df(arima_kalman_MAR_cor_50_errordf, ~as.data.frame(.x), .id="missingprop_autocor")

#60$
arima_kalman_MAR_cor_60_param<-purrr::map(arima_kalman_MAR_cor_60 , ~.["arima_pars"])
arima_kalman_MAR_cor_60_error<-purrr::map(arima_kalman_MAR_cor_60 , ~.["arima_errors"])

arima_kalman_MAR_cor_60_paramdf <- lapply(arima_kalman_MAR_cor_60_param, function(x) as.data.frame(do.call(rbind, x)))
arima_kalman_MAR_cor_60_errordf<- lapply(arima_kalman_MAR_cor_60_error, function(x) as.data.frame(do.call(rbind, x)))

kalman_MAR_cor_60_paramdf <- map_df(arima_kalman_MAR_cor_60_paramdf, ~as.data.frame(.x), .id="missingprop_autocor")
kalman_MAR_cor_60_errordf <- map_df(arima_kalman_MAR_cor_60_errordf, ~as.data.frame(.x), .id="missingprop_autocor")

#70$
arima_kalman_MAR_cor_70_param<-purrr::map(arima_kalman_MAR_cor_70 , ~.["arima_pars"])
arima_kalman_MAR_cor_70_error<-purrr::map(arima_kalman_MAR_cor_70 , ~.["arima_errors"])

arima_kalman_MAR_cor_70_paramdf <- lapply(arima_kalman_MAR_cor_70_param, function(x) as.data.frame(do.call(rbind, x)))
arima_kalman_MAR_cor_70_errordf<- lapply(arima_kalman_MAR_cor_70_error, function(x) as.data.frame(do.call(rbind, x)))

kalman_MAR_cor_70_paramdf <- map_df(arima_kalman_MAR_cor_70_paramdf, ~as.data.frame(.x), .id="missingprop_autocor")
kalman_MAR_cor_70_errordf <- map_df(arima_kalman_MAR_cor_70_errordf, ~as.data.frame(.x), .id="missingprop_autocor")

#80$
arima_kalman_MAR_cor_80_param<-purrr::map(arima_kalman_MAR_cor_80 , ~.["arima_pars"])
arima_kalman_MAR_cor_80_error<-purrr::map(arima_kalman_MAR_cor_80 , ~.["arima_errors"])

arima_kalman_MAR_cor_80_paramdf <- lapply(arima_kalman_MAR_cor_80_param, function(x) as.data.frame(do.call(rbind, x)))
arima_kalman_MAR_cor_80_errordf<- lapply(arima_kalman_MAR_cor_80_error, function(x) as.data.frame(do.call(rbind, x)))

kalman_MAR_cor_80_paramdf <- map_df(arima_kalman_MAR_cor_80_paramdf, ~as.data.frame(.x), .id="missingprop_autocor")
kalman_MAR_cor_80_errordf <- map_df(arima_kalman_MAR_cor_80_errordf, ~as.data.frame(.x), .id="missingprop_autocor")

#90#
arima_kalman_MAR_cor_90_param<-purrr::map(arima_kalman_MAR_cor_90 , ~.["arima_pars"])
arima_kalman_MAR_cor_90_error<-purrr::map(arima_kalman_MAR_cor_90 , ~.["arima_errors"])

arima_kalman_MAR_cor_90_paramdf <- lapply(arima_kalman_MAR_cor_90_param, function(x) as.data.frame(do.call(rbind, x)))
arima_kalman_MAR_cor_90_errordf<- lapply(arima_kalman_MAR_cor_90_error, function(x) as.data.frame(do.call(rbind, x)))

kalman_MAR_cor_90_paramdf <- map_df(arima_kalman_MAR_cor_90_paramdf, ~as.data.frame(.x), .id="missingprop_autocor")
kalman_MAR_cor_90_errordf <- map_df(arima_kalman_MAR_cor_90_errordf, ~as.data.frame(.x), .id="missingprop_autocor")

### rbind all the dataframes ## 

kalman_MAR_allcor_paramdf<-rbind(kalman_MAR_cor_01_paramdf, kalman_MAR_cor_10_paramdf, kalman_MAR_cor_20_paramdf, kalman_MAR_cor_30_paramdf, kalman_MAR_cor_40_paramdf, kalman_MAR_cor_50_paramdf, kalman_MAR_cor_60_paramdf,kalman_MAR_cor_70_paramdf,kalman_MAR_cor_80_paramdf,kalman_MAR_cor_90_paramdf)

rownames(kalman_MAR_allcor_paramdf) <- NULL

kalman_MAR_allcor_paramdf2<-kalman_MAR_allcor_paramdf %>% filter(!missingprop_autocor=="y_noMiss" ) %>% 
                                                      dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% 
                                                       select(missingprop_autocor, ar1, intercept, light, discharge)%>% 
                                                       mutate(type="Kalman filter") %>% 
                                                      separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16)

kalman_MAR_allcor_paramdf2$missingprop<-as.numeric(str_extract(kalman_MAR_allcor_paramdf2$missingprop1, "\\d+\\.*\\d*"))

kalman_MAR_allcor_paramdf2$autocorr<-as.numeric(str_extract(kalman_MAR_allcor_paramdf2$autoCorr1, "\\d+\\.*\\d*"))

# long form ##

paramkalmanlong <- kalman_MAR_allcor_paramdf2 %>% gather(param, value, ar1:discharge, factor_key=TRUE) %>% select(type, missingprop, autocorr, param, value)

## merge the real param values ##

paramkalmanlong2<-merge(paramkalmanlong, trueestdf)

paramkalmanlong3<-paramkalmanlong2 %>% mutate(abserror=abs(value-realvalue))


##make heat map! ##

mytheme<- theme_bw()+ theme(axis.line.x= element_line(colour = "black", size=0.3))+theme(axis.line.y= element_line(colour = "black", size=0.3))+theme(axis.text.x=element_text(size=10, colour = "black"))+theme(axis.text.y=element_text(size=10, colour = "black"))+theme(axis.title=element_text(size=12))+theme(plot.title=element_text(size=7) +theme(plot.title = element_text(hjust = 0.5)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(plot.title = element_text(margin=margin(0,0,5,0)))



####make all figures like this one####

heat_kalman<-ggplot(paramkalmanlong3, aes(x=missingprop, y=autocorr)) + 
  geom_tile(aes(fill=abserror), size=5) + 
  facet_wrap(~param)+
  scale_fill_viridis(begin=1, end=0)+
  xlab("Proportion missing")+
  ylab("Autocorrellation in missingness") +mytheme+ggtitle("GPP recover estimates with Kalman Filter")
