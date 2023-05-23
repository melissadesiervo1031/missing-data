
# Load packages
library(here)
library(tidyverse)
library(rstan)
library(brms)
library(Amelia)
library(forecast)
library(xts)
library(nlme)
library(lubridate)

# Reference A. Stears' code with helpful function for removing data
# makeMissing()
source("Functions/missing_data_functions.R")


#### Read in the missing dataframes that Alice S. made #####
###### Pull out the first in 1000 nested lists for this code ### (Eventually loop over all the lists)


gauss_sim_MAR_datasets <- readRDS(here("data/Missingdatasets/gauss_sim_randMiss.rds"))


##For nested list of GPP datasets with increasing MAR data add back in the covariates## 

GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]]

simmissingdf <-lapply(X = GPP_sim_MAR$y, 
                      FUN = function(X) cbind.data.frame(GPP = X, 
                                                         light = GPP_sim_MAR$sim_params$X[,2], 
                                                         discharge = GPP_sim_MAR$sim_params$X[,3]))



################################################################################
#### Data augmentation in STAN ###
################################################################################


# Create data list for feeding into Stan code:
# Generate vectors with indicies of observed and missing data points ##
# subset GPP timeseries to the observations only ##
stan_datasim <- lapply(simmissingdf, function(x) {
    ii_obs <- which(!is.na(x$GPP))
    ii_mis <- which(is.na(x$GPP))
    simdat <- list(
    N_obs = length(ii_obs),   # number of observations
    N_mis = length(ii_mis),   # number of missing data points
    ii_obs = ii_obs,          # indices of observations
    ii_mis = ii_mis,          # indices of missing observations
    P_obs = x$GPP[ii_obs],    # vector of observations
    light = x$light, discharge = x$discharge)
  
    return(simdat)
})

#### Model Fit ####

# Stan settings
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Fit models
# Fit datasets using BRMS:

bform <- brms::bf(GPP | mi() ~ light + discharge + ar(p = 1))
bmod <- brms::brm_multiple(bform, data = simmissingdf, combine = FALSE)

bfit <- bmod[[1]]
extract_brms_pars <- function(bfit){
    bsum <- summary(bfit$fit)$summary
    bsum <- data.frame(bsum) %>%
        mutate(parameter = row.names(bsum)) %>%
        select(parameter, mean, se_mean, sd, 
               '2.5%' = X2.5., '50%' = X50., '97.5%' = X97.5.) %>%
        slice(1:5)
  
    return(bsum)
}

bpars <- lapply(bmod, extract_brms_pars)
names(bpars) <- names(simmissingdf)
# can use this to see the estimated values of the missing data
# y_est <- brms::posterior_samples(bmod[[2]]$fit) %>%
#   summarize(across(everything(), mean)) %>%
#   select(starts_with('Ymi')) %>%
#   t() %>% as.vector()

# Fit models using Stan Missing data code
# Fit model
stan_datasim_fit <- lapply(stan_datasim,
                           function(x) stan("GPP sim and real/Stan_code/AR1_light_Q_centered_2.stan",
                                            data = x,
                                            chains = 4, 
                                            iter = 4000,
                                            control = list(max_treedepth = 12), 
                                            save_warmup=FALSE))



#stan_datasim_fit ## very large list ###


#### pull out SDs as well...###

##Pull param estimates into list ## takes a minute###
fit_summary_pars_bayes <- vector("list",16)
for (i in 1:16){
  fit_summary_pars_bayes[[i]]<-(summary(stan_datasim_fit[[i]], pars=c("beta[1]","beta[2]", "beta[3]", "phi","sigma"), probs=c(0.025,.5,.975))$summary)
}

names(fit_summary_pars_bayes) <- names(simmissingdf)

DAparamdf <- map_df(fit_summary_pars_bayes, ~as.data.frame(.x), .id="missingprop")
brmsparamdf <- map_df(bpars, ~as.data.frame(.x), .id="missingprop") %>%
  mutate(parameter = case_when(parameter == 'ar[1]'~ 'phi',
                               TRUE ~ parameter))

DAparamdf$parameter <- brmsparamdf$parameter

df <- left_join(DAparamdf, brmsparamdf, by = c('missingprop', 'parameter'), 
          suffix = c('', '.brms'))
true_pars <- data.frame(parameter = c('b_Intercept', 'b_light', 'b_discharge', 
                                      'phi', 'sigma'),
                        value = c(GPP_sim_MAR$sim_params$beta, 
                                  GPP_sim_MAR$sim_params$phi, 1))
df %>%
  mutate(missingprop = case_when(missingprop == 'y_noMiss' ~ 'propMissAct_0',
                                 TRUE ~ missingprop),
         missingprop = str_match(missingprop, 'propMissAct_([0-9.]+)$')[,2]) %>%
  left_join(true_pars, by = 'parameter') %>%
  mutate(diff.stan = mean - value, 
         diff.brms = mean.brms - value) %>%
  pivot_longer(cols = starts_with('diff'), 
               names_to = 'model', 
               values_to = 'diff', 
               names_pattern = 'diff.([a-z]+)') %>%
  ggplot(aes(missingprop, diff, col = model) )+
  geom_point() +
  facet_wrap(.~parameter) +
  theme_classic() +
  geom_hline(yintercept = 0)

# it looks like stan is still overestimating sigma and undersetimating phi...

############ MISSING NOT AT RANDOM MNAR ##############################

###### Pull out the first in 1000 nested lists for this code ### (Eventually loop over all the lists)

gauss_sim_MNAR_datasets <- readRDS(here("data/Missingdatasets/gauss_sim_minMaxMiss.rds"))

##For nested list of GPP datasets with increasing MNAR data add back in the date column and the covariates## 

GPP_sim_MNAR<- gauss_sim_MNAR_datasets [[1]]

simmissingNARdf <-lapply(X = GPP_sim_MNAR$y, 
                         FUN = function(X) cbind.data.frame(GPP = X, 
                                                            light = GPP_sim_MAR$sim_params$X[,2], 
                                                            discharge = GPP_sim_MAR$sim_params$X[,3]))

################################################################################
#### Data augmentation in STAN ###
################################################################################

# Create data list for feeding into Stan code:
# Generate vectors with indicies of observed and missing data points ##
# subset GPP timeseries to the observations only ##
stan_datasim_MNAR <- lapply(simmissingNARdf, function(x) {
    ii_obs <- which(!is.na(x$GPP))
    ii_mis <- which(is.na(x$GPP))
    simdat <- list(
    N_obs = length(ii_obs),   # number of observations
    N_mis = length(ii_mis),   # number of missing data points
    ii_obs = ii_obs,          # indices of observations
    ii_mis = ii_mis,          # indices of missing observations
    P_obs = x$GPP[ii_obs],    # vector of observations
    light = x$light, discharge = x$discharge)
  
  return(simdat)
})


#### Model Fit ####

# Stan settings
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Fit models
# Fit datasets using BRMS:

bform <- brms::bf(GPP | mi() ~ light + discharge + ar(p = 1))
bmodMNAR <- brms::brm_multiple(bform, data = simmissingNARdf, combine = FALSE)

bparsMNAR <- lapply(bmodMNAR, extract_brms_pars)
names(bparsMNAR) <- names(simmissingNARdf)
# can use this to see the estimated values of the missing data
# y_est <- brms::posterior_samples(bmod[[2]]$fit) %>%
#   summarize(across(everything(), mean)) %>%
#   select(starts_with('Ymi')) %>%
#   t() %>% as.vector()

# Fit models using Stan Missing data code
# Fit model
stan_datasim_MNAR_fit <- lapply(stan_datasim_MNAR,
                           function(x) stan("GPP sim and real/Stan_code/AR1_light_Q_centered_2.stan",
                                            data = x,
                                            chains = 4, 
                                            iter = 4000,
                                            control = list(max_treedepth = 12), 
                                            save_warmup=FALSE))


#### pull out SDs as well...###

##Pull param estimates into list ## takes a minute###
fit_summary_pars_bayesMNAR <- vector("list",16)
for (i in 1:16){
  fit_summary_pars_bayesMNAR[[i]]<-(summary(stan_datasim_MNAR_fit[[i]], 
                                        pars=c("beta[1]","beta[2]", "beta[3]", "phi","sigma"),
                                        probs=c(0.025,.5,.975))$summary)
}

names(fit_summary_pars_bayes_MNAR) <- names(simmissingNARdf)

DAparamdfMNAR <- map_df(fit_summary_pars_bayes_MNAR, ~as.data.frame(.x), .id="missingprop")
brmsparamdfMNAR <- map_df(bparsMNAR, ~as.data.frame(.x), .id="missingprop") %>%
  mutate(parameter = case_when(parameter == 'ar[1]'~ 'phi',
                               TRUE ~ parameter))

### formatting for figure ####

paramname<-c("intercept", "light", "discharge", "phi", "sdp")

param=rep(paramname, 20)

missingprop2=rep(missingprop, each=5)

DAparamdfMNAR2<-cbind(param=param, missingprop2=missingprop2, DAparamdfMNAR) 

DAparamdfMNAR3 <- DAparamdfMNAR2 %>% 
  mutate(type = "Data augmentation: STAN") %>% 
  select(missingprop2, type, param, mean, sd) %>% 
  rename(missingprop=missingprop2, value=mean, SD = sd) 
