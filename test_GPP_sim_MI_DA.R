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


## no missing###

sim1<-gauss_sim_MAR_datasets [[1]][["y"]][["y_noMiss"]]

realphi<-gauss_sim_MAR_datasets[[1]][["sim_params"]][["phi"]]
realbetas<-gauss_sim_MAR_datasets[[1]][["sim_params"]][["beta"]]

trueestdf <- data.frame (param = c("ar1", "intercept", "light", "discharge"), value = c(realphi, realbetas))

covariates<-gauss_sim_MAR_datasets[[1]][["sim_params"]][["X"]]

covariatesX<-as.matrix(covariates[,2:3])

days<-seq(1, 365)

sim1df<-as.data.frame(cbind(days=days, GPP=sim1, light=covariates[,2], discharge=covariates[,3]))


##quick plot simulated dataset##
simGPP <- ggplot(sim1df, aes(x=days, y=GPP))+
  geom_point(size = 2, color="chartreuse4") + 
  geom_line(size = 1, color="chartreuse4")+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title="Simulated GPP")



##For nested list of GPP datasets with increasing MAR data add back in the date column and the covariates## 

GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]][["y"]]

GPP_sim_MAR_2 <-lapply(X = GPP_sim_MAR, FUN = function(X)   cbind.data.frame(GPP=X, days=sim1df$days, light = sim1df$light, discharge = sim1df$discharge))



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


#### PLOT ALL THE ARIMA METHODS TOGETHER ###############


paramdroplong2 
paramMIlong2
paramNAlong2 

paramallARIMA<-rbind(paramdroplong2, paramMIlong2, paramNAlong2)


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

simmissingdf <- lapply(GPP_sim_MAR_2, function(x) as.data.frame(do.call(cbind, x)))

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

# Ran on server - started 3:27 , finished XXX.

#stan_datasim_fit ## very large list ###

##Pull param estimates into list ## takes a minute###
fit_summary_pars_bayes <- vector("list",20)
for (i in 1:20){
  fit_summary_pars_bayes[[i]]<-(summary(stan_datasim_fit[[i]], pars=c("beta[1]","beta[2]", "beta[3]", "phi","sdp"), probs=c(0.025,.5,.975))$summary)
}

names(fit_summary_pars_bayes) <- names(GPP_sim_MAR_2)


DAparamdf <- map_df(fit_summary_pars_bayes, ~as.data.frame(.x), .id="missingprop")

### formatting for figure ####

paramname<-c("intercept", "light", "discharge", "phi", "sdp")

param=rep(paramname, 20)

missingprop2=rep(missingprop, each=5)

DAparamdf2<-cbind(param=param, missingprop2=missingprop2, DAparamdf) 

DAparamdf3<- DAparamdf2 %>% mutate(type = "Data augmentation: STAN") %>% select(missingprop2, type, param, mean, se_mean) %>% rename(missingprop=missingprop2, value=mean, SE = se_mean) 



DApropmissingGPPsim<-ggplot(data=DAparamdf2, aes(x=as.numeric(missingprop2), y=mean))+
  facet_wrap(~param, scales="free",ncol=1)+
  geom_point(size=3)+
  #geom_hline(data=trueestdf, aes(yintercept=value), colour="salmon")+
  theme_classic()+
  ggtitle("Data augmentation: GPP simulated data")+
  xlab("Percent of Missing Data")+
  ylab("Parameter estimate")


#### plotting ARIMA AND DATA AUGMENTATION (STAN) TOGETHER ###

## change AR 1 to be phi ##

paramallARIMA2<-paramallARIMA

paramallARIMA2$param <- str_replace(paramallARIMA2$param, "ar1", "phi")

## distinguish between arima and stan ###

paramallARIMA2$type <- str_replace(paramallARIMA2$type, "Data Deletion", "Data Deletion: Arima")
paramallARIMA2$type <- str_replace(paramallARIMA2$type, "Multiple imputations", "Multiple imputations: Arima")
paramallARIMA2$type <- str_replace(paramallARIMA2$type, "Kalman filter", "Kalman filter: Arima (default)")

#### add ARIMA and STAN dataframes together ###

paramallARIMASTAN<-rbind(paramallARIMA2, DAparamdf3)

## drop sdp for comparison with arima ###
paramallARIMASTAN2<- paramallARIMASTAN %>% filter(param != "sdp" ) 

trueestdf2 <- data.frame (param = c("phi", "intercept", "light", "discharge"), value = c(realphi, realbetas))

arimastanMAR

paramallARIMASTAN2$type <- factor(paramallARIMASTAN2$type, levels=c("Data Deletion: Arima", "Kalman filter: Arima (default)", "Multiple imputations: Arima", "Data augmentation: STAN"), ordered=TRUE)

arimastanMAR<-ggplot(data=paramallARIMASTAN2, aes(x=as.numeric(missingprop), y=value))+
  facet_grid(~factor(param, levels=c("intercept", "phi", "light", "discharge"),exclude = NA)~ type, scales="free_y")+
  geom_hline(data=trueestdf2, aes(yintercept=value), colour="salmon")+
  #geom_hline(data=arimaestdf, aes(yintercept=value), colour="light blue")+
  geom_point(size=1)+
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), size=0.3)+
  theme_bw()+
  xlab("Percent of Missing Data (Missing at Random)")+
  ylab("Parameter estimate")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 6))
 
  
pdf(file = "Figures/MARcomparison.pdf",width = 8, height = 5)
arimastanMAR
dev.off()

######################## REPEAT FOR MNAR ###############################################################

###### Pull out the first in 1000 nested lists for this code ### (Eventually loop over all the lists)



gauss_sim_MNAR_datasets <- readRDS(here("data/Missingdatasets/gauss_sim_minMaxMiss.rds"))


## no missing###

sim1MNAR<-gauss_sim_MNAR_datasets [[1]][["y"]][["y_noMiss"]]

realphiMNAR<-gauss_sim_MNAR_datasets[[1]][["sim_params"]][["phi"]]
realbetasMNAR<-gauss_sim_MNAR_datasets[[1]][["sim_params"]][["beta"]]

trueestdfMNAR <- data.frame (param = c("ar1", "intercept", "light", "discharge"), value = c(realphiMNAR, realbetasMNAR))

covariates<-gauss_sim_MNAR_datasets[[1]][["sim_params"]][["X"]]

covariatesX<-as.matrix(covariates[,2:3])

days<-seq(1, 365)

sim1MNARdf<-as.data.frame(cbind(days=days, GPP=sim1MNAR, light=covariates[,2], discharge=covariates[,3]))

##For nested list of GPP datasets with increasing MNAR data add back in the date column and the covariates## 

GPP_sim_MNAR<- gauss_sim_MNAR_datasets [[1]][["y"]]

GPP_sim_MNAR_2 <-lapply(X = GPP_sim_MNAR, FUN = function(X)   cbind.data.frame(GPP=X, days=sim1df$days, light = sim1df$light, discharge = sim1df$discharge))



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

sim_missing_list_dropMNAR <- lapply(seq_along(GPP_sim_MNAR_2), function(j) {
  drop_na(GPP_sim_MNAR_2[[j]])
})

ArimaoutputdropMNAR <- lapply(seq_along(sim_missing_list_dropMNAR), function(j) {
  modeldrop <- Arima(sim_missing_list_dropMNAR[[j]][["GPP"]],order = c(1,0,0), xreg = matrix(c(sim_missing_list_dropMNAR[[j]][["light"]],sim_missing_list_dropMNAR[[j]][["discharge"]]), ncol = 2))
  arimacoefsdrop<-modeldrop$coef
  arimasesdrop<-sqrt(diag(vcov(modeldrop)))
  list(arimacoefsdrop=arimacoefsdrop, arimasesdrop=arimasesdrop)
})

names(ArimaoutputdropMNAR ) <- names(GPP_sim_MNAR_2)

########### formatting for figure #############

modeldropparamlistMNAR<-purrr::map(ArimaoutputdropMNAR , ~.["arimacoefsdrop"])
modeldropSElistMNAR<-map(ArimaoutputdropMNAR , ~.["arimasesdrop"])

modeldropparamlistMNAR2 <- lapply(modeldropparamlistMNAR, function(x) as.data.frame(do.call(rbind, x)))
modeldropSElistMNAR2 <- lapply(modeldropSElistMNAR, function(x) as.data.frame(do.call(rbind, x)))


modeldropparamdfMNAR <- map_df(modeldropparamlistMNAR2, ~as.data.frame(.x), .id="missingprop")
modeldropSEdfMNAR <- map_df(modeldropSElistMNAR2, ~as.data.frame(.x), .id="missingprop")

modeldropdfMNAR<-modeldropparamdfMNAR  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Data Deletion")

modeldropSEdfMNAR<-modeldropSEdfMNAR  %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Data Deletion")


missingprop<-seq(0.00, 0.95, by=0.05)

modeldropdfMNAR2<-cbind(missingprop, modeldropdfMNAR)

SEdropdfMNAR2<-cbind(missingprop, modeldropSEdfMNAR)


paramdroplongMNAR <- gather(modeldropdfMNAR2, param, value, ar1:discharge, factor_key=TRUE)

paramdropSElongMNAR <- gather(SEdropdfMNAR2, param, SE, ar1:discharge, factor_key=TRUE)

paramdroplongMNAR2<-merge(paramdroplongMNAR, paramdropSElongMNAR)

################################################## #########################################################
# ARIMA WITH NAS KALMAN FILTER
############################################################################################################

ArimaoutputNAsMNAR <- lapply(seq_along(GPP_sim_MNAR_2), function(j) {
  modelNAs <- Arima(GPP_sim_MNAR_2[[j]][["GPP"]],order = c(1,0,0), xreg = X)
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


################################################## #########################################################
# MULTIPLE IMPUTATIONS W/ AMELIA     ARIMA TO SOLVE FOR PARAMETER ESTIMATES       AVERAGE ESTIMATES
############################################################################################################

amelia1simMNAR <-lapply(X = GPP_sim_MNAR_2 , FUN = function(X)   amelia(X, ts="days", m=5, lags="GPP")) ## lags by 1 day ##


##nested list of dataframes that just has the imputations###
amelias11simMNAR<-map(amelia1simMNAR , ~.[["imputations"]])


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


#### PLOT ALL THE ARIMA METHODS TOGETHER ###############


paramdroplongMNAR2 
paramMIlongMNAR2
paramNAlongMNAR2 

paramallARIMAMNAR<-rbind(paramdroplongMNAR2, paramMIlongMNAR2, paramNAlongMNAR2)


arimallmethodsMNAR<-ggplot(data=paramallARIMAMNAR, aes(x=as.numeric(missingprop), y=value))+
  facet_grid(~factor(param, levels=c("ar1", "intercept", "light", "discharge"))~ type, scales="free")+
  geom_hline(data=trueestdf, aes(yintercept=value), colour="salmon")+
  geom_hline(data=arimaestdf, aes(yintercept=value), colour="light blue")+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE))+
  theme_bw()+
  ggtitle("Parameter estimates generated from Arima models")+
  xlab("Percent of Missing Data (Missing NOT at Random)")+
  ylab("Parameter estimate")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())





#####################################################################################################
#### Data augmentation in STAN ###
###################################################################################################

# Need to first decompose each of these nested lists into dfs

simmissingdfMNAR <- lapply(GPP_sim_MNAR_2, function(x) as.data.frame(do.call(cbind, x)))

# Also need to add sdo.
simmissingdfMNAR <- lapply(simmissingdfMNAR, function(x) cbind(x, sdo = 0.1))


# And a column to denote missingness and remove NAs from GPP data.
simmissingdfMNAR <- lapply(simmissingdfMNAR, function(x) x %>%
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

stan_datasimMNAR <- lapply(simmissingdfMNAR, function(x) stan_data_compile(x))

# Fit model
stan_datasim_fitMNAR <- lapply(stan_datasimMNAR,
                           function(x) stan("GPP sim and real/Stan_code/AR1_light_Q_centered.stan",
                                            data = x,
                                            chains = 4, 
                                            iter = 4000,
                                            control = list(max_treedepth = 12), 
                                            save_warmup=FALSE))

# Ran on server - started 5:10 , finished 5:38.

#stan_datasim_fit ## very large list ###

##Pull param estimates into list ## takes a minute###
fit_summary_pars_bayesMNAR <- vector("list",20)
for (i in 1:20){
  fit_summary_pars_bayesMNAR[[i]]<-(summary(stan_datasim_fitMNAR[[i]], pars=c("beta[1]","beta[2]", "beta[3]", "phi","sdp"), probs=c(0.025,.5,.975))$summary)
}

names(fit_summary_pars_bayesMNAR) <- names(GPP_sim_MNAR_2)


DAparamdfMNAR <- map_df(fit_summary_pars_bayesMNAR, ~as.data.frame(.x), .id="missingprop")

### formatting for figure ####

paramname<-c("intercept", "light", "discharge", "phi", "sdp")

param=rep(paramname, 20)

missingprop2=rep(missingprop, each=5)

DAparamdfMNAR2<-cbind(param=param, missingprop2=missingprop2, DAparamdfMNAR) 

DAparamdfMNAR3<- DAparamdfMNAR2 %>% mutate(type = "Data augmentation: STAN") %>% select(missingprop2, type, param, mean, se_mean) %>% rename(missingprop=missingprop2, value=mean, SE = se_mean) 


#### plotting ARIMA AND DATA AUGMENTATION (STAN) TOGETHER ###

## change AR 1 to be phi ##

paramallARIMAMNAR2<-paramallARIMAMNAR

paramallARIMAMNAR2$param <- str_replace(paramallARIMAMNAR2$param, "ar1", "phi")

## distinguish between arima and stan ###

paramallARIMAMNAR2$type <- str_replace(paramallARIMAMNAR2$type, "Data Deletion", "Data Deletion: Arima")
paramallARIMAMNAR2$type <- str_replace(paramallARIMAMNAR2$type, "Multiple imputations", "Multiple imputations: Arima")
paramallARIMAMNAR2$type <- str_replace(paramallARIMAMNAR2$type, "Kalman filter", "Kalman filter: Arima (default)")

#### add ARIMA and STAN dataframes together ###

paramallARIMASTANMNAR<-rbind(paramallARIMAMNAR2, DAparamdfMNAR3)

## drop sdp for comparison with arima ###
paramallARIMASTANMNAR2<- paramallARIMASTANMNAR %>% filter(param != "sdp" ) 

trueestdf2 <- data.frame (param = c("phi", "intercept", "light", "discharge"), value = c(realphi, realbetas))

arimastanMAR

paramallARIMASTANMNAR2$type <- factor(paramallARIMASTANMNAR2$type, levels=c("Data Deletion: Arima", "Kalman filter: Arima (default)", "Multiple imputations: Arima", "Data augmentation: STAN"), ordered=TRUE)

arimastanMNAR<-ggplot(data=paramallARIMASTANMNAR2, aes(x=as.numeric(missingprop), y=value))+
  facet_grid(~factor(param, levels=c("intercept", "phi", "light", "discharge"),exclude = NA)~ type, scales="free_y")+
  geom_hline(data=trueestdf2, aes(yintercept=value), colour="salmon")+
  #geom_hline(data=arimaestdf, aes(yintercept=value), colour="light blue")+
  geom_point(size=1)+
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), size=0.3)+
  theme_bw()+
  xlab("Percent of Missing Data (Missing NOT at Random)")+
  ylab("Parameter estimate")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 6))


pdf(file = "Figures/MNARcomparison.pdf",width = 8, height = 5)
arimastanMNAR
dev.off()


