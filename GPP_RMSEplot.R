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


############# MISSING AT RANDOM ########################## 

###pull 3 of Alice's nested lists of missing at random GPP data with diff param##

gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss.rds")

# nested list w/ 3 sims ##

GPP_3sims<- list(gauss_sim_MAR_datasets [[1]], gauss_sim_MAR_datasets [[2]], gauss_sim_MAR_datasets [[3]])

# add back in the date and X variables #

##For nested list of GPP datasets with increasing MAR data add back in the date column and the covariates## 

days<-seq(1, 365) 

covariateslist=purrr::map(GPP_3sims , ~.[["sim_params"]][["X"]])

lightlist=purrr::map(covariateslist, ~.[,2])

dischargelist=purrr::map(covariateslist, ~.[,3])


GPP_3sims2=list()
for (i in seq_along(GPP_3sims)) {
  a=list()
  for (j in seq_along(GPP_3sims[[i]][["y"]])) {
    tempobj=data.frame(days=days,GPP=GPP_3sims[[i]][["y"]][[j]], light=lightlist[i], discharge=dischargelist[i])
    colnames(tempobj)=c("days", "GPP", "light", "discharge")
    namelist <- names(GPP_3sims[[i]][["y"]])
    name=namelist[j]
    a[[name]] <- tempobj
  }
  
  GPP_3sims2[[i]] <- a
}



###


################### ARIMA method 1: Drop NAs ##################

source("Functions/Arima_drop_function.R")

### Run the ARIMA-DROP FUNCTION OVER THE 3 DATASETS###

arima_drop_MAR=list()
a=list()
for (i in seq_along(GPP_3sims)) {
  tempobj<- fit_arima_dropmissing(GPP_3sims[[i]]$y,GPP_3sims[[i]]$sim_params)
  names(tempobj)<-names(GPP_3sims[[i]][["y"]])
  arima_drop_MAR[[i]] <- tempobj
}

##name the elements in list for bookkeeping# 

names(arima_drop_MAR) <- paste("sim", 1:length(arima_drop_MAR), sep = "")

### pull out the output##

missingprop<-seq(0.00, 0.95, by=0.05)

modeldropparamlist<-list()
modeldropparamdf<-list()
paramdroplong<-list()

## estimated model parameters###

for (i in seq_along(arima_drop_MAR)) {
  
  tempobj<-purrr::map(arima_drop_MAR[[i]] , ~.["arima_pars"])
  modeldropparamlist[[i]]<-tempobj
  modeldropparamlist[[i]] <- lapply(modeldropparamlist[[i]], function(x) as.data.frame(do.call(rbind, x)))
  modeldropparamdf[[i]] <- map_df(modeldropparamlist[[i]], ~as.data.frame(.x), .id="missingprop")
  modeldropparamdf[[i]]<-modeldropparamdf[[i]]  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% mutate(missingprop=seq(0.00, 0.95, by=0.05)) %>% mutate(type="Data Deletion") %>% select(missingprop, ar1, intercept, light, discharge, type)
  rownames(modeldropparamdf[[i]]) <- NULL
  paramdroplong[[i]] <- gather(modeldropparamdf[[i]], param, value, ar1:discharge, factor_key=TRUE)
  
  }

##model errors###

modeldropSElist<-list()
modeldropSEdf<-list()
paramdropSElong<-list()

for (i in seq_along(arima_drop_MAR)) {
  
  tempobj2<-purrr::map(arima_drop_MAR[[i]] , ~.["arima_errors"])
  modeldropSElist[[i]]<-tempobj2
  modeldropSElist[[i]] <- lapply(modeldropSElist[[i]], function(x) as.data.frame(do.call(rbind, x)))
  modeldropSEdf[[i]] <- map_df(modeldropSElist[[i]], ~as.data.frame(.x), .id="missingprop")
  modeldropSEdf[[i]]<-modeldropSEdf[[i]]  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% mutate(missingprop=seq(0.00, 0.95, by=0.05)) %>% mutate(type="Data Deletion") %>% select(missingprop, ar1, intercept, light, discharge, type)
  rownames(modeldropSEdf[[i]])<-NULL
  paramdropSElong[[i]] <- gather(modeldropSEdf[[i]], param, SE, ar1:discharge, factor_key=TRUE)
}

#### pull out the REAL parameter values ####

paramlist<-list()

for (i in seq_along(arima_drop_MAR)) {
  tempobj1<-arima_drop_MAR[[i]][["y_noMiss"]][["sim_params"]][["phi"]]
  tempobj2<-arima_drop_MAR[[i]][["y_noMiss"]][["sim_params"]][["beta"]]
  tempobj3<-c(tempobj1, tempobj2)
  names<-c("ar1", "intercept", "light", "discharge")
  tempobj4<-cbind.data.frame(tempobj3, names)
  colnames(tempobj4)<-c("realvalue", "param")
  paramlist[[i]]<-tempobj4
}



## merge the parameters w/ SEs and the real values and calculate differences ###

paramdroplong2<-list()

for (i in seq_along(paramdroplong)) {
 tempobj3<-merge(paramdroplong[i], paramdropSElong[i])
 tempobj4<-merge(tempobj3, paramlist[i])
 paramdroplong2[[i]]<-tempobj4
 paramdroplong2[[i]]<-mutate(paramdroplong2[[i]], diff=realvalue-value, diffsquared=diff^2)
}


##name the elements in list for bookkeeping# 

names(paramdroplong2) <- paste("sim", 1:length(paramdroplong2), sep = "")

### Collapse the list and calculate RMSE ###

dropdf <- map_df(paramdroplong2, ~as.data.frame(.x), .id="id")

##summarize##

dropdfsummary<-dropdf%>% group_by(param, missingprop, type) %>% dplyr::summarise(sumdiffsquared=sum(diffsquared), n=n())%>% mutate(MSE=sumdiffsquared/(n-1))%>% mutate(RMSE=sqrt(MSE))


################## ARIMA method 2: KALMAN FILTER ##################

source("Functions/Arima_Kalman_function.R")

### Run the ARIMA-KALMAN FUNCTION OVER THE 3 DATASETS###

arima_kalman_MAR=list()
a=list()
for (i in seq_along(GPP_3sims)) {
  tempobj<- fit_arima_Kalman(GPP_3sims[[i]]$y,GPP_3sims[[i]]$sim_params)
  names(tempobj)<-names(GPP_3sims[[i]][["y"]])
  arima_kalman_MAR[[i]] <- tempobj
}


##name the elements in list for bookkeeping# 

names(arima_kalman_MAR) <- paste("sim", 1:length(arima_kalman_MAR), sep = "")

### pull out the output##

missingprop<-seq(0.00, 0.95, by=0.05)

modelkalmanparamlist<-list()
modelkalmanparamdf<-list()
paramkalmanlong<-list()

## estimated model parameters###


for (i in seq_along(arima_kalman_MAR)) {
  
  tempobj<-purrr::map(arima_kalman_MAR[[i]] , ~.["arima_pars"])
  modelkalmanparamlist[[i]]<-tempobj
  modelkalmanparamlist[[i]] <- lapply(modelkalmanparamlist[[i]], function(x) as.data.frame(do.call(rbind, x)))
  modelkalmanparamdf[[i]] <- map_df(modelkalmanparamlist[[i]], ~as.data.frame(.x), .id="missingprop")
  modelkalmanparamdf[[i]]<-modelkalmanparamdf[[i]]  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% mutate(missingprop=seq(0.00, 0.95, by=0.05)) %>% mutate(type="Kalman Filter") %>% select(missingprop, ar1, intercept, light, discharge, type)
  rownames(modelkalmanparamdf[[i]]) <- NULL
  paramkalmanlong[[i]] <- gather(modelkalmanparamdf[[i]], param, value, ar1:discharge, factor_key=TRUE)
  
}

##model errors###

modelkalmanSElist<-list()
modelkalmanSEdf<-list()
paramkalmanSElong<-list()

for (i in seq_along(arima_kalman_MAR)) {
  
  tempobj2<-purrr::map(arima_kalman_MAR[[i]] , ~.["arima_errors"])
  modelkalmanSElist[[i]]<-tempobj2
  modelkalmanSElist[[i]] <- lapply(modelkalmanSElist[[i]], function(x) as.data.frame(do.call(rbind, x)))
  modelkalmanSEdf[[i]] <- map_df(modelkalmanSElist[[i]], ~as.data.frame(.x), .id="missingprop")
  modelkalmanSEdf[[i]]<-modelkalmanSEdf[[i]]  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% mutate(missingprop=seq(0.00, 0.95, by=0.05)) %>% mutate(type="Kalman Filter") %>% select(missingprop, ar1, intercept, light, discharge, type)
  rownames(modelkalmanSEdf[[i]])<-NULL
  paramkalmanSElong[[i]] <- gather(modelkalmanSEdf[[i]], param, SE, ar1:discharge, factor_key=TRUE)
}

#### pull out the REAL parameter values ####  ## THIS SHOULD BE THE SAME FOR ALL METHODS..DOUBLE CHECK IF ITS NOT#

paramlist2<-list()

for (i in seq_along(arima_kalman_MAR)) {
  tempobj1<-arima_kalman_MAR[[i]][["y_noMiss"]][["sim_params"]][["phi"]]
  tempobj2<-arima_kalman_MAR[[i]][["y_noMiss"]][["sim_params"]][["beta"]]
  tempobj3<-c(tempobj1, tempobj2)
  names<-c("ar1", "intercept", "light", "discharge")
  tempobj4<-cbind.data.frame(tempobj3, names)
  colnames(tempobj4)<-c("realvalue", "param")
  paramlist2[[i]]<-tempobj4
}


## merge the parameters w/ SEs and the real values and calculate differences ###

paramkalmanlong2<-list()

for (i in seq_along(paramkalmanlong)) {
  tempobj3<-merge(paramkalmanlong[i], paramkalmanSElong[i])
  tempobj4<-merge(tempobj3, paramlist2[i])
  paramkalmanlong2[[i]]<-tempobj4
  paramkalmanlong2[[i]]<-mutate(paramkalmanlong2[[i]], diff=realvalue-value, diffsquared=diff^2)
}


##name the elements in list for bookkeeping# 

names(paramkalmanlong2) <- paste("sim", 1:length(paramkalmanlong2), sep = "")

### Collapse the list and calculate RMSE ###

kalmandf <- map_df(paramkalmanlong2, ~as.data.frame(.x), .id="id")

##summarize##

kalmandfsummary<-kalmandf%>% group_by(param, missingprop, type) %>% dplyr::summarise(sumdiffsquared=sum(diffsquared), n=n())%>% mutate(MSE=sumdiffsquared/(n-1))%>% mutate(RMSE=sqrt(MSE))



################## ARIMA method 3: MULTIPLE IMPUTATIONS ##################

source("Functions/Arima_MI_function.R")

### Run the ARIMA-MI FUNCTION OVER THE 3 DATASETS###

arima_MI_MAR=list()
a=list()
for (i in seq_along(GPP_3sims)) {
  tempobj<- fit_arima_MI(GPP_3sims[[i]][["y"]],GPP_3sims[[i]]$sim_params, imputationsnum = 5) #specify # imputations#
  arima_MI_MAR[[i]] <- tempobj
}

##stores the model parameters and SEs in seperate lists ##



## estimated model parameters###


modelMIparamlist<-list()
modelMIparamdf<-list()
paramMIlong<-list()

for (i in seq_along(arima_MI_MAR)) {
  
  tempobj<-purrr::pluck(arima_MI_MAR[[i]], 1)
  modelMIparamlist[[i]]<-tempobj
  modelMIparamlist[[i]] <- lapply(modelMIparamlist[[i]], function(x) as.data.frame(do.call(rbind, x)))
  modelMIparamlist[[i]]<-modelMIparamlist[[i]][- 1]  ### remove the list with no missing data###
  modelMIparamdf[[i]] <- map_df(modelMIparamlist[[i]], ~as.data.frame(.x), .id="missingprop")
  modelMIparamdf[[i]]<-modelMIparamdf[[i]]  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% mutate(missingprop=seq(0.05, 0.95, by=0.05)) %>% mutate(type="Multiple Imputations") %>% select(missingprop, ar1, intercept, light, discharge, type)
  rownames(modelMIparamdf[[i]]) <- NULL
  paramMIlong[[i]] <- gather(modelMIparamdf[[i]], param, value, ar1:discharge, factor_key=TRUE)
 
}


##model errors###

modelMISElist<-list()
modelMISEdf<-list()
paramMISElong<-list()

for (i in seq_along(arima_MI_MAR)) {
  
  tempobj2<-purrr::pluck(arima_MI_MAR[[i]], 2)
  modelMISElist[[i]]<-tempobj2
  modelMISElist[[i]] <- lapply(modelMISElist[[i]], function(x) as.data.frame(do.call(rbind, x)))
  modelMISElist[[i]]<-modelMISElist[[i]][- 1]  ### remove the list with no missing data###
  modelMISEdf[[i]] <- map_df(modelMISElist[[i]], ~as.data.frame(.x), .id="missingprop")
  modelMISEdf[[i]]<-modelMISEdf[[i]]  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% mutate(missingprop=seq(0.05, 0.95, by=0.05)) %>% mutate(type="Multiple Imputations") %>% select(missingprop, ar1, intercept, light, discharge, type)
  rownames(modelMISEdf[[i]])<-NULL
  paramMISElong[[i]] <- gather(modelMISEdf[[i]], param, SE, ar1:discharge, factor_key=TRUE)
}

#### pull out the REAL parameter values ####  ## THIS SHOULD BE THE SAME FOR ALL METHODS..DOUBLE CHECK IF ITS NOT#

### use the previous version, multiple imputations doesn't work without missing data##

## merge the parameters w/ SEs and the real values and calculate differences ###

paramMIlong2<-list()

for (i in seq_along(paramMIlong)) {
  tempobj3<-merge(paramMIlong[i], paramMISElong[i])
  tempobj4<-merge(tempobj3, paramlist2[i])
  paramMIlong2[[i]]<-tempobj4
  paramMIlong2[[i]]<-mutate(paramMIlong2[[i]], diff=realvalue-value, diffsquared=diff^2)
}


##name the elements in list for bookkeeping# 

names(paramMIlong2) <- paste("sim", 1:length(paramMIlong2), sep = "")

### Collapse the list and calculate RMSE ###

MIdf <- map_df(paramMIlong2, ~as.data.frame(.x), .id="id")

##summarize##

MIdfsummary<-MIdf%>% group_by(param, missingprop, type) %>% dplyr::summarise(sumdiffsquared=sum(diffsquared), n=n())%>% mutate(MSE=sumdiffsquared/(n-1))%>% mutate(RMSE=sqrt(MSE))


###### MERGE THE THREE ARIMA APPROACHES AND PLOT THE RMSE ####

RMSEallARIMA<-rbind(dropdfsummary, kalmandfsummary, MIdfsummary)


## change AR 1 to be phi ##

RMSEallARIMA<-paramallARIMA

RMSEallARIMA$param <- str_replace(RMSEallARIMA$param, "ar1", "phi")

RMSEarima<-ggplot(data=RMSEallARIMA, aes(x=as.numeric(missingprop), y=RMSE))+
  facet_grid(~factor(param, levels=c("intercept", "phi", "light", "discharge"),exclude = NA)~ type, scales="fixed")+
  geom_point(size=1)+
  theme_bw()+
  xlab("Percent of Missing Data (Missing at Random)")+
  ylab("Root mean square error (3 simulations)")






