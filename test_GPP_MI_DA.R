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

summary(arimafit)$coef

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



################## Multiple imputations w/ AMELIA ###################

amelia1 <-lapply(X = pine_missing_list_11 , FUN = function(X)   amelia(X, ts="date", m=5, polytime=1))

### Amelia makes 5 version of each imputed dataset for each item in the list ###


##pull out one list for forloop##

prop0.5list<-amelia1[["propMissIn_0.05; propMissAct_0.05"]][["imputations"]]


### ARIMA model to run on 1 imputed datasets###


###

X = matrix(c(pr$light.rel, pr$Q), ncol = 2)

modelOutput_list <- replicate(length(5), rep(NULL), simplify = FALSE)


for(i in 1:5){
  ## fit an ARIMA model to each AMELIA dataset ##
  arimaMI_i <-Arima(prop0.5list[[i]]$GPP, order = c(1,0,0), xreg = X)
  modelOutput_list[[i]] <- arimaMI_i$coef
}


###loop over all the lists ### Left to do...3/21/23##


#### how does ARIMA handle NAs ? ##

miss20<-pine_missing_list_11[["propMissIn_0.2; propMissAct_0.2"]]

X = matrix(c(pr$light.rel, pr$Q), ncol = 2)

arimafitfull <-  Arima(xts(pr$GPP, order.by=as.POSIXct(pr$date)), order = c(1,0,0), xreg = X)

arimafitmissing20<- arimafit <- Arima(xts(miss20$GPP, order.by=as.POSIXct(miss20$date)), order = c(1,0,0), xreg = X)


## unclear...I think it just skips over them? ##
#https://stats.stackexchange.com/questions/346225/fitting-arima-to-time-series-with-missing-values#
