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

brms_fit_MAR <- fit_brms_model(GPP_sim_MAR$y,
                               GPP_sim_MAR$sim_params,
                               include_missing = FALSE)


bpars<-brms_fit_MAR$brms_pars

brmsparamdf <- map_df(bpars, ~as.data.frame(.x), .id="missingprop_autocor") %>%
  mutate(parameter = case_when(parameter == 'ar[1]'~ 'phi', TRUE ~ parameter))  %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) 



## pull out missing prop and autocor values ##

brmsparamdf$missingprop<-as.numeric(str_extract(brmsparamdf$missingprop1, "\\d+\\.*\\d*"))

brmsparamdf$autocorr<-as.numeric(str_extract(brmsparamdf$autoCorr1, "\\d+\\.*\\d*"))


## fix names of params ##

brmsparamdf2<-brmsparamdf  %>% mutate(missingness="MAR low auto")%>% mutate(type="Data Augmentation: STAN")  %>% mutate(param=brmsparamdf$parameter) %>% mutate(value=brmsparamdf$mean)%>% mutate(SE=brmsparamdf$sd) %>% select(type, param, value, SE, missingprop, autocorr) 

brmsparamdf2$param <- str_replace(brmsparamdf2$param, "b_Intercept", "intercept")
brmsparamdf2$param <- str_replace(brmsparamdf2$param, "b_light", "light")
brmsparamdf2$param <- str_replace(brmsparamdf2$param, "b_discharge", "discharge")

#write.csv(brmsparamdf2, "C:/Users/Melissa/Dropbox/Academic stuff/MODELSCAPES/Missing data/BRMSMARlowauto.csv")

remove(brms_fit_MAR)

################################################## #########################################################
# MISSING AT RANDOM high autocor ##
############################################################################################################


gauss_sim_MAR_datasets_high <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_90.rds")
GPP_sim_MAR_high<- gauss_sim_MAR_datasets_high [[1]]

brms_fit_MAR_high <- fit_brms_model(GPP_sim_MAR_high$y,
                               GPP_sim_MAR_high$sim_params,
                               include_missing = FALSE)


bpars_high<-brms_fit_MAR_high$brms_pars

brmsparamdf_high <- map_df(bpars_high, ~as.data.frame(.x), .id="missingprop_autocor") %>%
  mutate(parameter = case_when(parameter == 'ar[1]'~ 'phi', TRUE ~ parameter))  %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) 



## pull out missing prop and autocor values ##

brmsparamdf_high$missingprop<-as.numeric(str_extract(brmsparamdf_high$missingprop1, "\\d+\\.*\\d*"))

brmsparamdf_high$autocorr<-as.numeric(str_extract(brmsparamdf_high$autoCorr1, "\\d+\\.*\\d*"))


## fix names of params ##

brmsparamdf_high2<-brmsparamdf_high %>% mutate(missingness="MAR high auto") %>% mutate(type="Data Augmentation: STAN")  %>% mutate(param=brmsparamdf_high$parameter) %>% mutate(value=brmsparamdf_high$mean)%>% mutate(SE=brmsparamdf_high$sd) %>% select(type, missingness,param, value, SE, missingprop, autocorr) 

brmsparamdf_high2$param <- str_replace(brmsparamdf_high2$param, "b_Intercept", "intercept")
brmsparamdf_high2$param <- str_replace(brmsparamdf_high2$param, "b_light", "light")
brmsparamdf_high2$param <- str_replace(brmsparamdf_high2$param, "b_discharge", "discharge")


#write.csv(brmsparamdf_high2, "C:/Users/Melissa/Dropbox/Academic stuff/MODELSCAPES/Missing data/BRMSMARhighauto.csv")

remove(brms_fit_MAR)

################################################## #########################################################
# MISSING NOT AT RANDOM MNAR
############################################################################################################


gauss_sim_MNAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_minMaxMiss.rds")
GPP_sim_MNAR<- gauss_sim_MAR_datasets [[1]]

brms_fit_MNAR <- fit_brms_model(GPP_sim_MNAR$y,
                               GPP_sim_MNAR$sim_params,
                               include_missing = FALSE)


bparsMNAR<-brms_fit_MNAR$brms_pars

brmsparamdf_MNAR <- map_df(bparsMNAR, ~as.data.frame(.x), .id="missingprop_autocor") %>%
  mutate(parameter = case_when(parameter == 'ar[1]'~ 'phi', TRUE ~ parameter))  %>%  separate(missingprop_autocor, into = c('missingprop1', 'autoCorr1'), sep = 16) 



## pull out missing prop and autocor values ##

brmsparamdf_MNAR$missingprop<-as.numeric(str_extract(brmsparamdf_MNAR$missingprop1, "\\d+\\.*\\d*"))

brmsparamdf_MNAR$autocorr<-as.numeric(str_extract(brmsparamdf_MNAR$autoCorr1, "\\d+\\.*\\d*"))


## fix names of params ##

brmsparamdf_MNAR2<-brmsparamdf_MNAR  %>% mutate(missingness="MNAR")%>% mutate(type="Data Augmentation: STAN")  %>% mutate(param=brmsparamdf_MNAR$parameter) %>% mutate(value=brmsparamdf_MNAR$mean)%>% mutate(SE=brmsparamdf_MNAR$sd) %>% select(type,missingness, param, value, SE, missingprop, autocorr) 

brmsparamdf_MNAR2$param <- str_replace(brmsparamdf_MNAR2$param, "b_Intercept", "intercept")
brmsparamdf_MNAR2$param <- str_replace(brmsparamdf_MNAR2$param, "b_light", "light")
brmsparamdf_MNAR2$param <- str_replace(brmsparamdf_MNAR2$param, "b_discharge", "discharge")

#write.csv(brmsparamdf_MNAR2, "C:/Users/Melissa/Dropbox/Academic stuff/MODELSCAPES/Missing data/BRMSMNAR.csv")

#BRMSMNAR<-read.csv("C:/Users/Melissa/Dropbox/Academic stuff/MODELSCAPES/Missing data/BRMSMNAR.csv")
#BRMSMARhighauto<-read.csv("C:/Users/Melissa/Dropbox/Academic stuff/MODELSCAPES/Missing data/BRMSMARhighauto.csv")
#BRMSMARlowauto<-read.csv("C:/Users/Melissa/Dropbox/Academic stuff/MODELSCAPES/Missing data/BRMSMARlowauto.csv")

#fix typo in low auto name #

#colnames(BRMSMARlowauto)<-c("X", "type", "missingness", "param", "value", "SE", "missingprop", "autocorr")


