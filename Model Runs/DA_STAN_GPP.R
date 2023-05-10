
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


##For nested list of GPP datasets with increasing MAR data add back in the date column and the covariates## 

GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]][["y"]]

GPP_sim_MAR_2 <-lapply(X = GPP_sim_MAR, FUN = function(X)   cbind.data.frame(GPP=X, days=sim1df$days, light = sim1df$light, discharge = sim1df$discharge))




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


#### pull out SDs as well...###

##Pull param estimates into list ## takes a minute###
fit_summary_pars_bayes <- vector("list",20)
for (i in 1:20){
  fit_summary_pars_bayes[[i]]<-(summary(stan_datasim_fit[[i]], pars=c("beta[1]","beta[2]", "beta[3]", "phi","sdp"), probs=c(0.025,.5,.975))$summary)
}

names(fit_summary_pars_bayes) <- names(GPP_sim_MAR_2)


DAparamdf <- map_df(fit_summary_pars_bayes, ~as.data.frame(.x), .id="missingprop")



############ MISSING NOT AT RANDOM MNAR ##############################

###### Pull out the first in 1000 nested lists for this code ### (Eventually loop over all the lists)

gauss_sim_MNAR_datasets <- readRDS(here("data/Missingdatasets/gauss_sim_minMaxMiss.rds"))

##For nested list of GPP datasets with increasing MNAR data add back in the date column and the covariates## 

GPP_sim_MNAR<- gauss_sim_MNAR_datasets [[1]][["y"]]

GPP_sim_MNAR_2 <-lapply(X = GPP_sim_MNAR, FUN = function(X)   cbind.data.frame(GPP=X, days=sim1df$days, light = sim1df$light, discharge = sim1df$discharge))


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

DAparamdfMNAR3<- DAparamdfMNAR2 %>% mutate(type = "Data augmentation: STAN") %>% select(missingprop2, type, param, mean, sd) %>% rename(missingprop=missingprop2, value=mean, SD = sd) 

