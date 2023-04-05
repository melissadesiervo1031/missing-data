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

##upload pine river dataset##
dat <- read_csv('data/NWIS_MissingTS_subset.csv')
mdat <- read_csv('data/NWIS_MissingTSinfo_subset.csv')

id <- mdat$site_name[4]
pr <- dat %>% filter(site_name == id) %>% select(date, GPP, light, Q, GPP.upper,GPP.lower) %>% mutate(Jdate= yday(date), light.rel = light/max(light))
pr<-as.data.frame(pr)

pr1<-pr %>% select(GPP)


##quick plot of pr GPP dataset##
pineriverGPP <- ggplot(pr, aes(date, GPP))+
    geom_point(size = 2, color="chartreuse4") + 
    geom_line(size = 1, color="chartreuse4")+
    labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title="Pine River GPP")

##### #  ARIMA and STAN model with NO missing data #### ###note that there actually are some missing dates..ignore for now##

### ARIMA ###

X = matrix(c(pr$light.rel, pr$Q), ncol = 2)

arimafit <- Arima(xts(pr$GPP, order.by=as.POSIXct(pr$date)), order = c(1,0,0), xreg = X)


se<-sqrt(diag(vcov(arimafit)))

arimafullcoefdf<-as.data.frame(t(as.data.frame(arimafit$coef)))

arimafullsedf<-as.data.frame(t(as.data.frame(sqrt(diag(vcov(arimafit))))))


#ar1:0.7139     intercept: 1.9943  light.rel= 5.2758  Q = -1.6564  #


### STAN #### 

# calculate observation error based on 95% CIs

sigma_obs <- (pr$GPP.upper - pr$GPP.lower)/3.92
pr <- pr %>%
  select(date, GPP, light, Q) %>%
  mutate(across(-date, ~zoo::na.approx(.)))%>%
  mutate(light.rel = light/max(light))


#Prep models ####
# model file: "model types/fixed_oi_light_centered.stan"
model_lq <- stan_model("GPP sim and real/Stan_code/AR1_light_Q_centered.stan")

#Create data object
data <- list(N = nrow(pr), P_obs = pr$GPP,
             sdo = sigma_obs, light = pr$light.rel, Q = pr$Q,
             miss_vec = rep(1, nrow(pr)))

#Run Stan
fit_lq <- rstan::sampling(object=model_lq, data = data,  
                          iter = 4000, chains = 4)

print(fit_lq,digits=5)

## beta1 (intercept) = 1.17345, beta2 = 2.93198 (light), beta 3 = -1.26293 (Q), phi = 0.62793, sdp = 1.17685 ###


###comparing ARIMA and STAN ### similar estimates for everything except light. SDP in STAN only ###



####################Apply missing data function to Pine River  GPP ############

pine_missing_list<-makeMissing(timeSeries = pr$GPP, typeMissing = "random")#  makes a list of GPP w/ increasing missingness## Missing at random from 0.5 - 95%

      ##add back in the date column and the covariates## 

pine_missing_list_11 <-lapply(X = pine_missing_list, FUN = function(X)   cbind.data.frame(GPP=X, date=as.Date(pr$date, format="%Y-%m-%d"), Q = pr$Q, light.rel=pr$light.rel))


###################### Data augmentation in STAN #####################

## turn lists into dataframes ##

pinemissingdfs <- lapply(pine_missing_list_11, function(x) as.data.frame(do.call(cbind, x)))

##add SDO ##

sdo<-(pr$GPP.upper - pr$GPP.lower)/3.92

pinemissingdfs_2 <- lapply(pinemissingdfs, function(x) cbind(x, sdo))

# And a column to denote missingness and remove NAs from GPP data.
pinemissingdfs_3 <- lapply(pinemissingdfs_2, function(x) x %>%
                    mutate(miss_vec = case_when(is.na(GPP) == TRUE ~ 0,
                                                TRUE ~ 1)) %>%
                    mutate(GPP_noNA = case_when(is.na(GPP) == TRUE ~ 0,
                                                TRUE ~ GPP)))
# Stan settings
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Compile data
stan_data_compile <- function(x){
  data <- list(N = length(x$GPP), # number of records
               P_obs = x$GPP_noNA, # simulated GPP w/o NAs
               light = x$light, # relativized light
               Q = x$Q,   # relativized discharge
               sdo = x$sdo,  # standard deviation of GPP estimates
               miss_vec = x$miss_vec) # vector of missingness
  return(data)
}

stan_data <- lapply(pinemissingdfs_3 , function(x) stan_data_compile(x))

# Fit model
GPPstanmissing <- lapply(stan_data,
                   function(x) stan("GPP sim and real/Stan_code/AR1_light_Q_centered.stan",
                                    data = x,
                                    chains = 4, 
                                    iter = 4000,
                                    control = list(max_treedepth = 12), 
                                    save_warmup=FALSE))




# Ran on my laptop - started 5:02, finished 5:22.

saveRDS(GPPstanmissing, file="GPPstanmissing.RData")

outputlist<-print(GPPstanmissing, pars=c("phi", "sdp", "beta"), digits=7)

## need to pull items out of list to make figure##



################## Multiple imputations w/ AMELIA ARIMA models ###################

amelia1 <-lapply(X = pine_missing_list_11 , FUN = function(X)   amelia(X, ts="date", m=5, polytime=1))

### Amelia makes 5 version of each imputed dataset for each item in the list ###

##nested list of dataframes that just has the imputations###
amelias11<-map(amelia1 , ~.[["imputations"]])


###loop over all the lists ##

##start with amelia11, nested list ###

##matrix of covariates##
X = matrix(c(pr$light.rel, pr$Q), ncol = 2)

##working forloop## gives us the model parameters and errors for all the ARIMA models on imputed datasets

modelparamlist=list()
modelerrorlist=list()
for (i in seq_along(amelias11)) {
  a=list()
  aa=list()
  for (j in seq_along(amelias11[[i]])) {
    tempobj=Arima(amelias11[[i]][[j]]$GPP, order = c(1,0,0), xreg = X)
    arimacoefs<-tempobj$coef
    arimases<-sqrt(diag(vcov(tempobj)))
    name <- paste('imp',seq_along((amelias11)[[i]])[[j]],sep='')
    a[[name]] <- arimacoefs
    aa[[name]]<-arimases
  }
  name1 <- names(amelias11)[[i]]
  modelparamlist[[name1]] <- a
  modelerrorlist[[name1]] <- aa
}

modelparamlist
modelerrorlist


### Averages the models together back to 1 model per missing data prop ##

listcoefses<-mapply(function(X,Y) {
         list(mi.meld(data.frame(X), data.frame(Y), byrow=FALSE))
             }, X=modelparamlist, Y=modelerrorlist)

##pulls out parameters and ses ##

modelavgparamlist<-map(listcoefses , ~.["q.mi"])
modelseparamlist<-map(listcoefses , ~.["se.mi"])


modelavgparamdf <- map_df(modelavgparamlist, ~as.data.frame(.x), .id="missingprop")
modelSEdf <- map_df(modelseparamlist, ~as.data.frame(.x), .id="missingprop")

missingprop<-seq(from=0.05, to =0.95, by=0.05)

modelavgparamdf2<-modelavgparamdf %>% dplyr::rename(ar1=q.mi.ar1, intercept=q.mi.intercept, light=q.mi.xreg1, discharge=q.mi.xreg2) %>% select(ar1, intercept, light, discharge) %>% mutate(type="Multiple imputations")

modelavgparamdf2<-cbind(missingprop, modelavgparamdf2)

      ###add in the no missing data##

arimafullcoefdf2<-arimafullcoefdf %>% dplyr::rename(ar1=ar1, intercept=intercept, light=xreg1, discharge=xreg2) %>%   mutate(type="Multiple imputations")  %>% mutate(missingprop=0.00) %>% select(propmissing, ar1, intercept, light, discharge, type)

modelavgparamdf3<-rbind(as.data.frame(arimafullcoefdf2), as.data.frame(modelavgparamdf2))

######## RUN ARIMA with DELETED data, no imputation or Augmentation ##

head(pine_missing_list_11)


ArimaoutputNAs <- lapply(seq_along(pine_missing_list_11), function(j) {
  modelNAs <- Arima(pine_missing_list_11[[j]][["GPP"]],order = c(1,0,0), xreg = X)
  arimacoefsNAs<-modelNAs$coef
  arimasesNAs<-sqrt(diag(vcov(modelNAs)))
  list(arimacoefsNAs=arimacoefsNAs, arimasesNAs=arimasesNAs)
 })


names(ArimaoutputNAs) <- names(pine_missing_list_11)


modelNAparamlist<-map(ArimaoutputNAs , ~.["arimacoefsNAs"])
modelNASElist<-map(ArimaoutputNAs , ~.["arimasesNAs"])

modelNAparamlist2 <- lapply(modelNAparamlist, function(x) as.data.frame(do.call(rbind, x)))
modelNASElist2 <- lapply(modelNASElist, function(x) as.data.frame(do.call(rbind, x)))


modelNAparamdf <- map_df(modelNAparamlist2, ~as.data.frame(.x), .id="missingprop")
modelNASEdf <- map_df(modelNASElist2, ~as.data.frame(.x), .id="missingprop")

missingprop<-seq(from=0.05, to =0.95, by=0.05)



##



