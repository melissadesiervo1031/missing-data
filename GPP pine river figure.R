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
source("functions/missing_data_functions.R")

########## Run code for ARIMA methods over GPP 1 simulated dataset w/ increasing missingness
## 3 different types (MAR low auto, MAR high auto, and MNAR) #####

## ARIMA method 1: Drop NAs ###

source("Functions/Arima_drop_function.R")


## ARIMA method 2: Kalman filter (preserves the NAS) ###

source("Functions/Arima_Kalman_function.R")


## ARIMA method 3: Multiple imputations with AMELIA ###


source("Functions/Arima_MI_function.R")


########## Run code for ARIMA methods over GPP 1 simulated dataset w/ increasing missingness #####

source("Functions/model_fitting.R")


################################################## #########################################################
# MISSING AT RANDOM low autocor ##
############################################################################################################

PineRiverlowauto <- readRDS("data/Missingdatasets/gauss_real_randMiss_autoCorr_10.rds")

PineRiverlowauto_list<-as.list(PineRiverlowauto[,9:23])

Arimanomissing_PineRiver=arima(PineRiverlowauto$GPP, order = c(1,0,0), xreg = X)

PineRivercoef<-Arimanomissing_PineRiver[["coef"]]

phi=PineRivercoef[1]

betas=PineRivercoef[2:4]

X = matrix(c(rep(1, times=length(PineRiverlowauto$date)), PineRiverlowauto$light.rel,  PineRiverlowauto$Q), ncol = 3)

sim_pars_PR = list(phi=phi, betas=betas, X=X)

#drop #

arima_drop_pineriver_lowauto<- fit_arima_dropmissing(PineRiverlowauto_list, sim_pars=sim_pars_PR)

#kalman#

Kalman_pineriver_lowauto<- fit_arima_Kalman(PineRiverlowauto_list, sim_pars=sim_pars_PR)

#MI#

MI_pineriver_lowauto<- fit_arima_MI(PineRiverlowauto_list, sim_pars=sim_pars_PR)

# not working #



##ran them seperately and saved .csv. combining them all together back in here #

