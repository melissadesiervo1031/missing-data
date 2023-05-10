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
# LIST WISE DELETION SOLVE WITH ARIMA
############################################################################################################

# drops rows with missing data ###

sim_missing_list_drop <- lapply(seq_along(GPP_sim_MAR_2), function(j) {
  drop_na(GPP_sim_MAR_2[[j]])
})

Arimaoutputdrop <- lapply(seq_along(sim_missing_list_drop ), function(j) {
  modeldrop <- Arima(sim_missing_list_drop [[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(sim_missing_list_drop [[j]][["light"]],sim_missing_list_drop [[j]][["discharge"]]), ncol = 2))
  arimacoefsdrop<-modeldrop$coef
  arimasesdrop<-sqrt(diag(vcov(modeldrop)))
  list(arimacoefsdrop=arimacoefsdrop, arimasesdrop=arimasesdrop)
})

names(Arimaoutputdrop ) <- names(GPP_sim_MAR_2)

########### formatting for figure #############

modeldropparamlist<-purrr::map(Arimaoutputdrop , ~.["arimacoefsdrop"])
modeldropSElist<-map(Arimaoutputdrop , ~.["arimasesdrop"])

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