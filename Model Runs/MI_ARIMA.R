# Load packages
library(here)
library(tidyverse)
library(Amelia)
library(forecast)
library(xts)
library(nlme)
library(tidyverse)
library(lubridate)


# Reference A. Stears' code with helpful function for removing data
# makeMissing()
source("Functions/missing_data_functions.R")


#### Read in the missing dataframes that Alice S. made #####
###### Pull out the first in 1000 nested lists for this code ### (Eventually loop over all the lists)


###### Pull out the first in 1000 nested lists for this code ### (Eventually loop over all the lists)


gauss_sim_MAR_datasets <- readRDS(here("data/Missingdatasets/gauss_sim_randMiss.rds"))


##For nested list of GPP datasets with increasing MAR data add back in the date column and the covariates## 


GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]][["y"]]

sim1<-gauss_sim_MAR_datasets [[1]][["y"]][["y_noMiss"]]

covariates<-gauss_sim_MAR_datasets[[1]][["sim_params"]][["X"]]

covariatesX<-as.matrix(covariates[,2:3])

days<-seq(1, 365)


sim1df<-as.data.frame(cbind(days=days, GPP=sim1, light=covariates[,2], discharge=covariates[,3]))

GPP_sim_MAR_2 <-lapply(X = GPP_sim_MAR, FUN = function(X)   cbind.data.frame(GPP=X, days=sim1df$days, light = sim1df$light, discharge = sim1df$discharge))


#### MISSING COMPLETELY AT RANDOM (MCAR) ######

################################################## #########################################################
# MULTIPLE IMPUTATIONS W/ AMELIA     ARIMA TO SOLVE FOR PARAMETER ESTIMATES       AVERAGE ESTIMATES
############################################################################################################

amelia1sim <-lapply(X = GPP_sim_MAR_2 , FUN = function(X)   amelia(X, ts="days", m=5, lags="GPP")) ## lags by 1 day ##


##nested list of dataframes that just has the imputations###
amelias11sim<-map(amelia1sim , ~.[["imputations"]])


###
##matrix of covariates##
X = matrix(c(sim1df$light,  sim1df$discharge), ncol = 2)

##forloop## gives us the model parameters and errors for all the ARIMA models on imputed datasets

modelparamlistsim=list()
modelerrorlistsim=list()
for (i in seq_along(amelias11sim)) {
  a=list()
  aa=list()
  for (j in seq_along(amelias11sim[[i]])) {
    tempobj=Arima(amelias11sim[[i]][[j]]$GPP, order = c(1,0,0), xreg = X)
    arimacoefs<-tempobj$coef
    arimases<-sqrt(diag(vcov(tempobj)))
    name <- paste('imp',seq_along((amelias11sim)[[i]])[[j]],sep='')
    a[[name]] <- arimacoefs
    aa[[name]]<-arimases
  }
  name1 <- names(amelias11sim)[[i]]
  modelparamlistsim[[name1]] <- a
  modelerrorlistsim[[name1]] <- aa
}

modelparamlistsim
modelerrorlistsim


### Averages the models together back to 1 model per missing data prop ##

listcoefsessim<-mapply(function(X,Y) {
  list(mi.meld(data.frame(X), data.frame(Y), byrow=FALSE))
}, X=modelparamlistsim, Y=modelerrorlistsim)

##pulls out parameters and ses ##

paramlistsim<-map(listcoefsessim , ~.["q.mi"])
selistsim<-map(listcoefsessim , ~.["se.mi"])


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


###### Pull out the first in 1000 nested lists for this code ### (Eventually loop over all the lists)

gauss_sim_MNAR_datasets <- readRDS(here("data/Missingdatasets/gauss_sim_minMaxMiss.rds"))

##For nested list of GPP datasets with increasing MNAR data add back in the date column and the covariates## 

GPP_sim_MNAR<- gauss_sim_MNAR_datasets [[1]][["y"]]

GPP_sim_MNAR_2 <-lapply(X = GPP_sim_MNAR, FUN = function(X)   cbind.data.frame(GPP=X, days=sim1df$days, light = sim1df$light, discharge = sim1df$discharge))


################################################## #########################################################
# MULTIPLE IMPUTATIONS W/ AMELIA     ARIMA TO SOLVE FOR PARAMETER ESTIMATES       AVERAGE ESTIMATES
############################################################################################################

amelia1simMNAR <-lapply(X = GPP_sim_MNAR_2 , FUN = function(X)   amelia(X, ts="days", m=5, lags="GPP")) ## lags by 1 day ##


##nested list of dataframes that just has the imputations###
amelias11simMNAR<-map(amelia1simMNAR , ~.[["imputations"]])



###### Pull out the first in 1000 nested lists for this code ### (Eventually loop over all the lists)

gauss_sim_MNAR_datasets <- readRDS(here("data/Missingdatasets/gauss_sim_minMaxMiss.rds"))

## no missing###

sim1MNAR<-gauss_sim_MNAR_datasets [[1]][["y"]][["y_noMiss"]]

covariates<-gauss_sim_MNAR_datasets[[1]][["sim_params"]][["X"]]

covariatesX<-as.matrix(covariates[,2:3])

days<-seq(1, 365)

sim1dfMNAR<-as.data.frame(cbind(days=days, GPP=sim1MNAR, light=covariates[,2], discharge=covariates[,3]))



###
##matrix of covariates##
X = matrix(c(sim1dfMNAR$light,  sim1dfMNAR$discharge), ncol = 2)

##forloop## gives us the model parameters and errors for all the ARIMA models on imputed datasets

modelparamlistsimMNAR=list()
modelerrorlistsimMNAR=list()
for (i in seq_along(amelias11simMNAR)) {
  a=list()
  aa=list()
  for (j in seq_along(amelias11simMNAR[[i]])) {
    tempobj=Arima(amelias11simMNAR[[i]][[j]]$GPP, order = c(1,0,0), xreg = X)
    arimacoefs<-tempobj$coef
    arimases<-sqrt(diag(vcov(tempobj)))
    name <- paste('imp',seq_along((amelias11simMNAR)[[i]])[[j]],sep='')
    a[[name]] <- arimacoefs
    aa[[name]]<-arimases
  }
  name1 <- names(amelias11simMNAR)[[i]]
  modelparamlistsimMNAR[[name1]] <- a
  modelerrorlistsimMNAR[[name1]] <- aa
}

modelparamlistsimMNAR
modelerrorlistsimMNAR


### Averages the models together back to 1 model per missing data prop ##

listcoefsessimMNAR<-mapply(function(X,Y) {
  list(mi.meld(data.frame(X), data.frame(Y), byrow=FALSE))
}, X=modelparamlistsimMNAR, Y=modelerrorlistsimMNAR)

##pulls out parameters and ses ##

paramlistsimMNAR<-map(listcoefsessimMNAR , ~.["q.mi"])
selistsimMNAR<-map(listcoefsessimMNAR , ~.["se.mi"])


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
