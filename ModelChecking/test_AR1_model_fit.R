# Test Script for AR1 River GPP Model
# Created: March 15, 2023
# Edited by: Heili Lowman

# The following script is designed to be a test run to fit
# a few of the datasets created by the "gaussian_ar1_data_sims.R"
# script (stored locally in the "data" folder as "gauss_ar1_0miss_
# datasets.rds) using the AR1 model found in the "GPP sim and real/
# Stan_code" folder as "AR1_light_Q_centered.stan".

#### Setup ####

# Load packages
library(here)
library(tidyverse)
library(rstan)

# Reference A. Stears' code with helpful function for removing data
# makeMissing()
source("missing_data_functions.R")

# Load raw data & trim to site of interest.
raw_dat <- read_csv('data/NWIS_MissingTS_subset.csv')
sable_dat <- raw_dat %>%
  filter(site_name == "nwis_04137500")

# Load simulated data & trim to 5 trial time series.
dat <- readRDS("data/gauss_ar1_0miss_datasets.rds")
dat5 <- dat[c(4,38,650,843,999)]

# Quick plots of timeseries as they are
plot(x = 1:365, y = dat5[[1]]$y)
plot(x = 1:365, y = dat5[[2]]$y)
plot(x = 1:365, y = dat5[[3]]$y)
plot(x = 1:365, y = dat5[[4]]$y)
plot(x = 1:365, y = dat5[[5]]$y)

# Going to remove data using the "random chunks" approach since this
# most closely mimics removal of data following storm events (~14 days).

# Note - only accepts "X" as input variable
dat5_rc <- lapply(X = dat5, 
                  # function will apply to each of the 5 timeseries
                  FUN = function(X) makeMissing(timeSeries = X$y, 
                  # and will remove data in chunks
                                                typeMissing = "randChunks",
                  # proportion missing
                  # default is a vector from .05:.95 by .05    
                                                propMiss = 0.35,
                  # integer length of data missing in each chunk
                  # default is 5% of length of full time series
                                                chunkSize = 14))

# Quick plots of timeseries with 35% missingness in random chunks
# So, there should be roughly 128 days missing in 2 week-long chunks
# so, 9 chunks missing from each dataset.
plot(x = 1:365, y = dat5_rc[[1]]$propMissing_0.35)
plot(x = 1:365, y = dat5_rc[[2]]$propMissing_0.35)
plot(x = 1:365, y = dat5_rc[[3]]$propMissing_0.35)
plot(x = 1:365, y = dat5_rc[[4]]$propMissing_0.35)
plot(x = 1:365, y = dat5_rc[[5]]$propMissing_0.35)

# Looks good! Some chunks are put together while others aren't.

# Now, to add light and discharge to each of these GPP estimates.
sable_dat <- sable_dat %>%
  mutate(light.rel = light/max(light),
         Q.rel = scale(Q)[,1])

Q.rel <- sable_dat$Q.rel
light.rel <- sable_dat$light.rel

# Need to first decompose each of these nested lists into dfs
# because the joins are giving me issues.
dat5_rc <- lapply(dat5_rc, function(x) as.data.frame(do.call(cbind, x)))

# Add discharge & light data.
dat5_rc <- lapply(dat5_rc, function(x) cbind(x, Q = Q.rel))
dat5_rc <- lapply(dat5_rc, function(x) cbind(x, light = light.rel))

# Also need to add sdo.
dat5_rc <- lapply(dat5_rc, function(x) cbind(x, sdo = 0.1))

# And a column to denote missingness and remove NAs from GPP data.
dat5_rc <- lapply(dat5_rc, function(x) x %>%
    mutate(miss_vec = case_when(is.na(propMissing_0.35) == TRUE ~ 0,
                                                TRUE ~ 1)) %>%
      mutate(GPP_noNA = case_when(is.na(propMissing_0.35) == TRUE ~ 0,
                                  TRUE ~ propMissing_0.35)))

#### Model Fit ####

# Stan settings
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Compile data
stan_data_compile <- function(x){
  data <- list(N = length(x$propMissing_0.35), # number of records
               P_obs = x$GPP_noNA, # simulated GPP w/o NAs
               light = x$light, # relativized light
               Q = x$Q,   # relativized discharge
               sdo = x$sdo,  # standard deviation of GPP estimates
               miss_vec = x$miss_vec) # vector of missingness
  return(data)
}

stan_data5 <- lapply(dat5_rc, function(x) stan_data_compile(x))

# Fit model
ar1_fit5 <- lapply(stan_data5,
                   function(x) stan("GPP sim and real/Stan_code/AR1_light_Q_centered.stan",
                                    data = x,
                                    chains = 4, 
                                    iter = 4000,
                     control = list(max_treedepth = 12), 
                     save_warmup=FALSE))

# Ran on server - started 9:38, finished 9:42.

# Examine output
traceplot(ar1_fit5[[1]], pars=c("phi", "sdp", "beta"))
pairs(ar1_fit5[[1]], pars=c("phi", "sdp","beta","lp__"))
plot(ar1_fit5[[1]], pars=c("phi", "sdp", "beta"))
print(ar1_fit5[[1]], pars=c("phi", "sdp", "beta"))

# Trace plots and Rhat values look alright, and no divergent
# transitions spat out as warnings.

#         mean se_mean   sd  2.5%   25%  50%  75% 97.5% n_eff Rhat
# phi     0.95    0.00 0.02  0.91  0.94 0.95 0.97  0.99  3031    1
# sdp     0.86    0.00 0.04  0.78  0.84 0.86 0.89  0.95  4898    1
# beta[1] 2.88    0.03 2.28 -2.34  1.67 3.03 4.25  6.95  4371    1
# beta[2] 0.04    0.00 0.16 -0.28 -0.07 0.04 0.14  0.34  7462    1
# beta[3] 0.04    0.00 0.05 -0.06  0.00 0.04 0.07  0.14  6342    1

# End of script.