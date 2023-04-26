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



gauss_miss_datasets <- readRDS(here("data/Missingdatasets/gauss_sim_randMiss.rds"))



#########################################################################################################################
## Pull in simulated AR time series and create lists of datasets w/ NAS (missing at random)
###########################################################################################################################


##Pull in one simulated GPP dataset to try out methods##

gauss_ar1_0miss_datasets <- readRDS(here("data/gauss_ar1_0miss_datasets.rds"))

sim1<-gauss_ar1_0miss_datasets[[1]]$y 

#### params #####
realphi<-gauss_ar1_0miss_datasets[[1]]$sim_params$phi
realbetas<-gauss_ar1_0miss_datasets[[1]]$sim_params$beta

trueestdf <- data.frame (param = c("ar1", "intercept", "light", "discharge"), value = c(realphi, realbetas))


####

covariates<-gauss_ar1_0miss_datasets[[1]]$sim_params$X

covariatesX<-as.matrix(covariates[,2:3])

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

#########################################################################################################################
## ARIMA estimate with no missing data ##
###########################################################################################################################

##matrix of covariates##
X = matrix(c(sim1df$light,  sim1df$discharge), ncol = 2)

Arimanomissing=arima(sim1df$GPP, order = c(1,0,0), xreg = X)

Arimanomissing

arimaestdf <- data.frame (param = c("ar1", "intercept", "light", "discharge"), value =Arimanomissing$coef)


################################################## #########################################################
# LIST WISE DELETION SOLVE WITH ARIMA
############################################################################################################

# drops rows with missing data ###

sim_missing_list_drop <- lapply(seq_along(sim_missing_list_2), function(j) {
 drop_na(sim_missing_list_2[[j]])
})

Arimaoutputdrop <- lapply(seq_along(sim_missing_list_drop ), function(j) {
  modeldrop <- Arima(sim_missing_list_drop [[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(sim_missing_list_drop [[j]][["light"]],sim_missing_list_drop [[j]][["discharge"]]), ncol = 2))
  arimacoefsdrop<-modeldrop$coef
  arimasesdrop<-sqrt(diag(vcov(modeldrop)))
  list(arimacoefsdrop=arimacoefsdrop, arimasesdrop=arimasesdrop)
})

names(Arimaoutputdrop ) <- names(sim_missing_list_2)



modeldropparamlist<-purrr::map(Arimaoutputdrop , ~.["arimacoefsdrop"])
modeldropSElist<-map(Arimaoutputdrop , ~.["arimasesdrop"])

modeldropparamlist2 <- lapply(modeldropparamlist, function(x) as.data.frame(do.call(rbind, x)))
modeldropSElist2 <- lapply(modeldropSElist, function(x) as.data.frame(do.call(rbind, x)))


modeldropparamdf <- map_df(modeldropparamlist2, ~as.data.frame(.x), .id="missingprop")
modeldropSEdf <- map_df(modeldropSElist2, ~as.data.frame(.x), .id="missingprop")

modeldropdf<-modeldropparamdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Data Deletion")

modeldropSEdf<-modeldropSEdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Data Deletion")


missingprop<-seq(0.05, 0.95, by=0.05)

modeldropdf2<-cbind(missingprop, modeldropdf)

SEdropdf2<-cbind(missingprop, modeldropSEdf)


paramdroplong <- gather(modeldropdf2, param, value, ar1:discharge, factor_key=TRUE)

paramdropSElong <- gather(SEdropdf2, param, SE, ar1:discharge, factor_key=TRUE)

paramdroplong2<-merge(paramdroplong, paramdropSElong)

deletepropmissingGPPsim<-ggplot(data=paramdroplong, aes(x=as.numeric(missingprop), y=value))+
  facet_wrap(~param, scales="free",ncol=1)+
  geom_point(size=3)+
  geom_hline(data=trueestdf, aes(yintercept=value), colour="salmon")+
  theme_classic()+
  ggtitle("Data deletion: GPP simulated data")+
  xlab("Percent of Missing Data")+
  ylab("Parameter estimate")

################################################## #########################################################
# ARIMA WITH NAS KALMAN FILTER
############################################################################################################

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

modelNAdf<-modelNAparamdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Kalman filter")

modelNASEdf<-modelNASEdf  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Kalman filter")


missingprop<-seq(0.05, 0.95, by=0.05)

modelNAdf2<-cbind(missingprop, modelNAdf)

SENAdf2<-cbind(missingprop, modelNASEdf)


paramNAlong <- gather(modelNAdf2, param, value, ar1:discharge, factor_key=TRUE)

paramNASElong <- gather(SENAdf2, param, SE, ar1:discharge, factor_key=TRUE)

paramNAlong2<-merge(paramNAlong, paramNASElong)


################################################## #########################################################
# MULTIPLE IMPUTATIONS W/ AMELIA     ARIMA TO SOLVE FOR PARAMETER ESTIMATES       AVERAGE ESTIMATES
############################################################################################################

amelia1sim <-lapply(X = sim_missing_list_2 , FUN = function(X)   amelia(X, ts="days", m=5, lags="GPP")) ## lags by 1 day ##

###plot the time series with missing ##

simmissing40<-sim_missing_list_2[[8]]


##quick plot simulated dataset with missing ##
plotmissing40 <- ggplot(simmissing40, aes(x=days, y=GPP))+
  geom_point(size = 2, color="chartreuse4") + 
  geom_line(size = 1, color="chartreuse4")+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title="Simulated GPP with 40 % missing data")

########## plot what amelia is imputing ####

a.out8<-amelia1sim[8] 

a.outimp1<-a.out8$`propMissIn_0.4; propMissAct_0.4`$imputations$imp1
missmatrix<-a.out8$`propMissIn_0.4; propMissAct_0.4`$missMatrix[,1]

imp1missing40<-cbind(a.outimp1, miss=missmatrix)

##quick plot what amelia is imputing ##
plotimpmissing40 <- ggplot(imp1missing40, aes(x=days, y=GPP, color=miss))+
  geom_point(size = 2) + 
  geom_line(size = 1, color="chartreuse4")+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title="Simulated GPP with 40 % missing data, imputed values")



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

avglSEdf2<-avglSEdf  %>% dplyr::rename(ar1=se.mi.ar1, intercept=se.mi.intercept, light=se.mi.xreg1, discharge=se.mi.xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Multiple imputations")


missingprop<-seq(0.05, 0.95, by=0.05)

MIdf2<-cbind(missingprop, avgparamdf2)

MISEdf2<-cbind(missingprop, avglSEdf2)


paramMIlong <- gather(MIdf2, param, value, ar1:discharge, factor_key=TRUE)

paramMISElong <- gather(MISEdf2, param, SE, ar1:discharge, factor_key=TRUE)

paramMIlong2<-merge(paramMIlong,paramMISElong)


head(Arimanomissingdf)


#### PLOT ALL THE ARIMA METHODS TOGETHER ###############


paramdroplong2 
paramMIlong2
paramNAlong2 

paramallARIMA<-rbind(paramdroplong2, paramMIlong2, paramNAlong2)

cust_label <- setNames(paste0("sch.id:", unique(ten$sch.id)), unique(ten$sch.id))

arimallmethods<-ggplot(data=paramallARIMA, aes(x=as.numeric(missingprop), y=value))+
  facet_grid(~factor(param, levels=c("ar1", "intercept", "light", "discharge"))~ type, scales="free")+
  geom_hline(data=trueestdf, aes(yintercept=value), colour="salmon")+
  geom_hline(data=arimaestdf, aes(yintercept=value), colour="light blue")+
    geom_point(size=3)+
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE))+
  theme_bw()+
  ggtitle("Parameter estimates generated from Arima models")+
  xlab("Percent of Missing Data (Missing at Random)")+
  ylab("Parameter estimate")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())





#####################################################################################################
#### Data augmentation in STAN ###
###################################################################################################

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

##Pull param estimates into list
fit_summary_pars_bayes <- vector("list",19)
for (i in 1:19){
  fit_summary_pars_bayes[[i]]<-(summary(stan_datasim_fit[[i]], pars=c("beta[1]","beta[2]", "beta[3]", "phi","sdp"), probs=c(0.025,.5,.975))$summary)
}

names(fit_summary_pars_bayes) <- names(sim_missing_list_2)


DAparamdf <- map_df(fit_summary_pars_bayes, ~as.data.frame(.x), .id="missingprop")

paramname<-c("beta[1]", "beta[2]", "beta[3]", "phi", "sdp")

param=rep(paramname, 19)

missingprop2=rep(missingprop, each=5)

DAparamdf2<-cbind(param=param, missingprop2=missingprop2, DAparamdf) 


DApropmissingGPPsim<-ggplot(data=DAparamdf2, aes(x=as.numeric(missingprop2), y=mean))+
  facet_wrap(~param, scales="free",ncol=1)+
  geom_point(size=3)+
  #geom_hline(data=trueestdf, aes(yintercept=value), colour="salmon")+
  theme_classic()+
  ggtitle("Data augmentation: GPP simulated data")+
  xlab("Percent of Missing Data")+
  ylab("Parameter estimate")



