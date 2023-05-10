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


# Reference A. Stears' code with helpful function for removing data
# makeMissing()
source("missing_data_functions.R")


#### Read in the missing dataframes that Alice S. made #####
###### Pull out the first in 1000 nested lists for this code ### (Eventually loop over all the lists)


gauss_sim_MAR_datasets <- readRDS(here("data/Missingdatasets/gauss_sim_randMiss.rds"))


#### MISSING COMPLETELY AT RANDOM (MCAR) ######


################################################## #########################################################
# ARIMA WITH NAS KALMAN FILTER
############################################################################################################

ArimaoutputNAs <- lapply(seq_along(GPP_sim_MAR_2), function(j) {
  modelNAs <- Arima(GPP_sim_MAR_2[[j]][["GPP"]],order = c(1,0,0), xreg = X)
  arimacoefsNAs<-modelNAs$coef
  arimasesNAs<-sqrt(diag(vcov(modelNAs)))
  list(arimacoefsNAs=arimacoefsNAs, arimasesNAs=arimasesNAs)
})

names(ArimaoutputNAs) <- names(GPP_sim_MAR_2)

############ formatting for figure #############

modelNAparamlist<-purrr::map(ArimaoutputNAs , ~.["arimacoefsNAs"])
modelNASElist<-map(ArimaoutputNAs , ~.["arimasesNAs"])

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
