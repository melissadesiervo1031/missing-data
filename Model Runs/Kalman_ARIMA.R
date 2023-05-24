# Load packages
library(here)
library(tidyverse)
library(Amelia)
library(forecast)
library(xts)
library(nlme)
library(tidyverse)
library(lubridate)


#### Read in the missing dataframes that Alice S. made #####
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
# ARIMA WITH NAS KALMAN FILTER
############################################################################################################

ArimaoutputNAs <- lapply(seq_along(GPP_sim_MAR_2), function(j) {
  modelNAs <- Arima(GPP_sim_MAR_2[[j]][["GPP"]],order = c(1,0,0), xreg = covariatesX)
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



############ MISSING NOT AT RANDOM MNAR ##############################


###### Pull out the first in 1000 nested lists for this code ### (Eventually loop over all the lists)

gauss_sim_MNAR_datasets <- readRDS(here("data/Missingdatasets/gauss_sim_minMaxMiss.rds"))

##For nested list of GPP datasets with increasing MNAR data add back in the date column and the covariates## 

GPP_sim_MNAR<- gauss_sim_MNAR_datasets [[1]][["y"]]

GPP_sim_MNAR_2 <-lapply(X = GPP_sim_MNAR, FUN = function(X)   cbind.data.frame(GPP=X, days=sim1df$days, light = sim1df$light, discharge = sim1df$discharge))


################################################## #########################################################
# ARIMA WITH NAS KALMAN FILTER
############################################################################################################

ArimaoutputNAsMNAR <- lapply(seq_along(GPP_sim_MNAR_2), function(j) {
  modelNAs <- Arima(GPP_sim_MNAR_2[[j]][["GPP"]],order = c(1,0,0), xreg = covariatesX)
  arimacoefsNAs<-modelNAs$coef
  arimasesNAs<-sqrt(diag(vcov(modelNAs)))
  list(arimacoefsNAs=arimacoefsNAs, arimasesNAs=arimasesNAs)
})

names(ArimaoutputNAsMNAR) <- names(GPP_sim_MNAR_2)

############ formatting for figure #############

modelNAparamlistMNAR<-purrr::map(ArimaoutputNAsMNAR , ~.["arimacoefsNAs"])
modelNASElistMNAR<-map(ArimaoutputNAsMNAR , ~.["arimasesNAs"])

modelNAparamlistMNAR2 <- lapply(modelNAparamlistMNAR, function(x) as.data.frame(do.call(rbind, x)))
modelNASElistMNAR2 <- lapply(modelNASElistMNAR, function(x) as.data.frame(do.call(rbind, x)))


modelNAparamdfMNAR <- map_df(modelNAparamlistMNAR2, ~as.data.frame(.x), .id="missingprop")
modelNASEdfMNAR <- map_df(modelNASElistMNAR2, ~as.data.frame(.x), .id="missingprop")

modelNAdfMNAR<-modelNAparamdfMNAR  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Kalman filter")

modelNASEdfMNAR<-modelNASEdfMNAR  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Kalman filter")


missingprop<-seq(0.00, 0.95, by=0.05)

modelNAdfMNAR2<-cbind(missingprop, modelNAdfMNAR)

SENAdfMNAR2<-cbind(missingprop, modelNASEdfMNAR)


paramNAlongMNAR <- gather(modelNAdfMNAR2, param, value, ar1:discharge, factor_key=TRUE)

paramNASElongMNAR <- gather(SENAdfMNAR2, param, SE, ar1:discharge, factor_key=TRUE)

paramNAlongMNAR2<-merge(paramNAlongMNAR, paramNASElongMNAR)
