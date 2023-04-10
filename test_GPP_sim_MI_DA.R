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

##Pull in one simulated GPP dataset to try out methods##

gauss_ar1_0miss_datasets <- readRDS("C:/Users/Melissa/Dropbox/Academic stuff/MODELSCAPES/Missing data/missing-data-git/data/gauss_ar1_0miss_datasets.rds")

sim1<-gauss_ar1_0miss_datasets[[1]]$y 

covariates<-gauss_ar1_0miss_datasets[[1]]$sim_params$X

days<-seq(1, 365)

sim1df<-as.data.frame(cbind(days=days, GPP=sim1, light=covariates[,2], discharge=covariates[,3]))


##quick plot simulated dataset##
simGPP <- ggplot(sim1df, aes(x=days, y=GPP))+
  geom_point(size = 2, color="chartreuse4") + 
  geom_line(size = 1, color="chartreuse4")+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title="Simulated GPP")


###make lists of missing data##

sim_missing_list<-makeMissing(timeSeries = sim1df$GPP, propMiss <- seq(0.05, 0.95, by = .05), typeMissing = "random")#  makes a list of GPP w/ increasing missingness## Missing at random from 0.5 - 95%

##add back in the date column and the covariates## 

sim_missing_list_2 <-lapply(X = sim_missing_list, FUN = function(X)   cbind.data.frame(GPP=X, days=sim1df$days, light = sim1df$light, discharge = sim1df$discharge))


## ARIMA estimate with no missing data ##

Arimanomissing=Arima(sim_missing_list_2$`propMissIn_0; propMissAct_0`$GPP, order = c(1,0,0), xreg = X)

Arimanomissing

Arimanomissingdf<-as.data.frame(t(as.data.frame(Arimanomissing$coef)))



#### data augmentation in STAN ###

# Need to first decompose each of these nested lists into dfs

simmissingdf <- lapply(sim_missing_list_2, function(x) as.data.frame(do.call(cbind, x)))

# Also need to add sdo.
simmissingdf <- lapply(simmissingdf, function(x) cbind(x, sdo = 0.1))


# And a column to denote missingness and remove NAs from GPP data.
simmissingdf <- lapply(simmissingdf, function(x) x %>%
                    mutate(miss_vec = case_when(is.na(GPP) == TRUE ~ 0,
                                                TRUE ~ 1)) %>%
                    mutate(GPP_noNA = case_when(is.na(GPP) == TRUE ~ 0,
                                                TRUE ~ GPP)))

#### Model Fit ####

# Stan settings
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Compile data
stan_data_compile <- function(x){
  data <- list(N = length(x$GPP), # number of records
               P_obs = x$GPP_noNA, # simulated GPP w/o NAs
               light = x$light, # relativized light
               Q = x$discharge,   # relativized discharge
               sdo = x$sdo,  # standard deviation of GPP estimates
               miss_vec = x$miss_vec) # vector of missingness
  return(data)
}

stan_datasim <- lapply(simmissingdf , function(x) stan_data_compile(x))

# Fit model
stan_datasim_fit <- lapply(stan_datasim,
                   function(x) stan("GPP sim and real/Stan_code/AR1_light_Q_centered.stan",
                                    data = x,
                                    chains = 4, 
                                    iter = 4000,
                                    control = list(max_treedepth = 12), 
                                    save_warmup=FALSE))

# Ran on server - started 12:52 , finished XXX.

#stan_datasim_fit ## very large list ###

paramDA<-map(stan_datasim_fit , ~.["samples"])



################## Multiple imputations w/ AMELIA ARIMA models ###################

amelia1sim <-lapply(X = sim_missing_list_2 , FUN = function(X)   amelia(X, ts="days", m=5, polytime=1))

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

missingprop<-seq(from=0.05, to =0.95, by=0.05)

avgparamdf2<-avgparamdf %>% dplyr::rename(ar1=q.mi.ar1, intercept=q.mi.intercept, light=q.mi.xreg1, discharge=q.mi.xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Multiple imputations")

avgparamdf2<-cbind(missingprop, avgparamdf2)

paramMIlong <- gather(avgparamdf2, param, value, ar1:discharge, factor_key=TRUE)

### calculate error in parameter estimates ##

realar1=gauss_ar1_0miss_datasets[[1]]$sim_params$phi
realint=2.6806 ###i took the value from the arima intercept..not sure we need it##
realB1=-1.17030279 ###from the simulation#
realB2=0.05385505 ###from the simulation#


avgparamdf2sim<-avgparamdf2 %>% mutate(ar1pererror= (ar1-realar1)/realar1) %>% mutate(intpererror= (intercept-realint)/realint)%>% mutate(lightpererror= (light-realB1)/realB1)%>% mutate(dischargepererror= (discharge-realB2)/realB2)%>% mutate(meanerror=round(rowMeans(abs(avgparamdf2sim[,7:10]), na.rm=TRUE),digits=3))

      ####GRAPH THE ACTUAL PARAMETER ESTIMATES###

head(Arimanomissingdf)

trueestdf <- data.frame (param = c("ar1", "intercept", "light", "discharge"), value = c(0.9798821,2.680645,-1.112911,0.04300555))


  mipropmissingGPPsim<-ggplot(data=paramMIlong, aes(x=as.numeric(missingprop), y=value))+
  facet_wrap(~param, scales="free",ncol=1)+
  geom_point(size=3)+
  geom_hline(data=trueestdf, aes(yintercept=value), colour="salmon")+
  theme_classic()+
  ggtitle("Multiple Imputations: GPP simulated data")+
  xlab("Percent of Missing Data")+
  ylab("Parameter estimate")


#### DELETING DATA ###
  
  ArimaoutputNAs <- lapply(seq_along(sim_missing_list_2), function(j) {
    modelNAs <- Arima(sim_missing_list_2[[j]][["GPP"]],order = c(1,0,0), xreg = X)
    arimacoefsNAs<-modelNAs$coef
    arimasesNAs<-sqrt(diag(vcov(modelNAs)))
    list(arimacoefsNAs=arimacoefsNAs, arimasesNAs=arimasesNAs)
  })

  names(ArimaoutputNAs) <- names(sim_missing_list_2)
  
  
  modelNAparamlist<-purrr::map(ArimaoutputNAs , ~.["arimacoefsNAs"])
  modelNASElist<-map(ArimaoutputNAs , ~.["arimasesNAs"])
  
  modelNAparamlist2 <- lapply(modelNAparamlist, function(x) as.data.frame(do.call(rbind, x)))
  modelNASElist2 <- lapply(modelNASElist, function(x) as.data.frame(do.call(rbind, x)))
  

  modelNAparamdf <- map_df(modelNAparamlist2, ~as.data.frame(.x), .id="missingprop")
  modelNASEdf <- map_df(modelNASElist2, ~as.data.frame(.x), .id="missingprop")
  
  modelNAdf<-modelNAparamdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Data deletion")
  
  modelNAdf2<-cbind(missingprop, modelNAdf)
  
  paramNAlong <- gather(modelNAdf2, param, value, ar1:discharge, factor_key=TRUE)
  
  deletepropmissingGPPsim<-ggplot(data=paramNAlong, aes(x=as.numeric(missingprop), y=value))+
    facet_wrap(~param, scales="free",ncol=1)+
    geom_point(size=3)+
    geom_hline(data=trueestdf, aes(yintercept=value), colour="salmon")+
    theme_classic()+
    ggtitle("Data Deletion: GPP simulated data")+
    xlab("Percent of Missing Data")+
    ylab("Parameter estimate")
  