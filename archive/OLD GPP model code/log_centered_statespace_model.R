# setwd('C:/Users/alice.carter/git/missing-data/')

library(rstan)
library(lubridate)
library(tidyverse)
library(streamMetabolizer)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#Simulate data

##light
# Use the built in stream Metabolizer function:
# this function gives instantaneous light throughout the day. For our purposes,
# we want daily light. To get this, we will grab light at 15 minute intervals 
# then sum over each day.

# calculate at Flathead Bio Station:
lat = 47.8762
lon = -114.03

light_inst <- data.frame(datetime_MST = seq(from = ymd_hms('2020-01-01 00:00:00', tz = 'MST'),
                                            to = ymd_hms('2020-12-31 23:45:00', tz = 'MST'), 
                                            by = '15 min')) %>%
  mutate(solartime = streamMetabolizer::calc_solar_time(datetime_MST, lon),
         light = streamMetabolizer::calc_light(solartime, lat, lon))

# calculate daily light 

light_daily <- light_inst %>%
  mutate(date = as.Date(datetime_MST, tz = 'MST')) %>%
  group_by(date) %>%
  summarize(light = sum(light)) %>%
  ungroup() %>%
  mutate(light.rel = light/max(light))

light.rel = light_daily$light.rel

N <- nrow(light_daily) #length of data
z <- numeric(N + 1)
tt <- 1:N

##GPP based on time-series model with known parameters
# GPP_t = b0 + phi*GPP_t-1 + b1*L + proc_error
# GPP_obs = GPP + obs_error
b0 <- 3      # intercept
b1 <- 1   # light coefficient
phi <- 0.5   # ar1 coefficient
sigma_proc <- 0.3 # ~5% error
sigma_obs <- 0.1  # fixed in model

# Expected value based on parameters:  ~10.4
GPP <- vector(mode = 'double', length = N)  # GPP state variable
GPP[1] <- b0 + b1 * light.rel[1] + 
                rnorm(1, 0, sd = sigma_proc/sqrt(1-phi^2))

# Observed GPP: add observation error on a linear, not a log scale
P_obs = rnorm(1, GPP[1], sigma_obs)
for (i in 2:N){
  GPP[i] <- b0 + b1 * light.rel[i] + phi * (GPP[i-1] - b0) + 
                  rnorm(1, 0, sigma_proc)
  P_obs[i] = rnorm(1, GPP[i], sigma_obs)
}

plot(GPP, type = 'l', ylab = 'GPP')
points(P_obs)
legend('topright', c('state', 'observations'), pch = c(NA, 1), lty = c(1, 0),
       bty = 'n', ncol = 2)

##Bind simulated data and time
y_full <- as.data.frame(cbind(GPP, P_obs, doy = as.factor(tt)))


# model file: "model types/fixed_oi_light_centered.stan"
#Prep model
model <- "model types/fixed_oi_light_centered.stan"
model <- stan_model(model)

#Create data object
data <- list(N = length(P_obs), P_obs = P_obs,
             sdo = sigma_obs, light=light.rel)

#Run Stan
fit <- rstan::sampling(object=model, data = data,  
                       iter = 4000, chains = 4)

##Print and extract
fit_extract <- rstan::extract(fit)
print(fit, pars=c("phi", "sdp", "beta" ))
pairs(fit, pars=c("phi", "sdp","beta","lp__"))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit)

##Traceplots
traceplot(fit, pars=c("phi", "sdp", "beta"))

##Plot density plots
plot_sdp <- stan_dens(fit, pars="sdp") + geom_vline(xintercept = sigma_proc)
plot_phi <- stan_dens(fit, pars="phi") + geom_vline(xintercept = phi)
plot_b0 <- stan_dens(fit, pars="beta[1]") + geom_vline(xintercept = b0) +
  xlab("Intercept")
plot_b1 <- stan_dens(fit, pars="beta[2]") + geom_vline(xintercept = b1) +
  xlab("Light beta")
gridExtra::grid.arrange(plot_b0, plot_b1, plot_phi, plot_sdp, nrow=2)

# simulate data for model with light and discharge
# get some real discharge data:
qq <- read_csv('data/NWIS_MissingTS_subset.csv')
qq <- filter(qq, site_name == qq$site_name[1]) %>%
  bind_rows(qq) %>%
  slice(1:N)

##GPP based on time-series model with known parameters
# GPP_t = b0 + phi*GPP_t-1 + b1*L + proc_error
# GPP_obs = GPP + obs_error
b0 <- 3      # intercept
b1 <- 1   # light coefficient
b2 <- -1   # Q coefficient
phi <- 0.5   # ar1 coefficient
sigma_proc <- 0.3 # ~5% error
sigma_obs <- 0.1  # fixed in model

# Expected value based on parameters:  ~10.4
GPP <- vector(mode = 'double', length = N)  # GPP state variable
GPP[1] <- b0 + b1 * light.rel[1] + b2*qq$Q[1] +
                rnorm(1, 0, sd = sigma_proc/sqrt(1-phi^2))

# Observed GPP: add observation error on a linear, not a log scale
P_obs = rnorm(1, GPP[1], sigma_obs)
for (i in 2:N){
  GPP[i] <- b0 + b1 * light.rel[i] + b2 * qq$Q[i] + phi * (GPP[i-1] - b0) + 
                  rnorm(1, 0, sigma_proc)
  P_obs[i] = rnorm(1, GPP[i], sigma_obs)
}

plot(GPP, type = 'l', ylab = 'GPP')
points(P_obs)
legend('topright', c('state', 'observations'), pch = c(NA, 1), lty = c(1, 0),
       bty = 'n', ncol = 2)

##Bind simulated data and time
y_full <- as.data.frame(cbind(GPP, P_obs, doy = as.factor(tt)))


# model file: "model types/fixed_oi_light_Q_centered.stan"
#Prep model
model <- "model types/fixed_oi_light_Q_centered.stan"
model <- stan_model(model)

#Create data object
data <- list(N = length(P_obs), P_obs = P_obs, Q = qq$Q,
             sdo = sigma_obs, light=light.rel)

#Run Stan
fit <- rstan::sampling(object=model, data = data,  
                       iter = 4000, chains = 4)

##Print and extract
print(fit, pars=c("phi", "sdp", "beta" ))
pairs(fit, pars=c("phi", "sdp","beta","lp__"))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit)

##Traceplots
traceplot(fit, pars=c("phi", "sdp", "beta"))

##Plot density plots
plot_sdp <- stan_dens(fit, pars="sdp") + geom_vline(xintercept = sigma_proc)
plot_phi <- stan_dens(fit, pars="phi") + geom_vline(xintercept = phi)
plot_b0 <- stan_dens(fit, pars="beta[1]") + geom_vline(xintercept = b0) +
  xlab("Intercept")
plot_b1 <- stan_dens(fit, pars="beta[2]") + geom_vline(xintercept = b1) +
  xlab("Light beta")
plot_b2 <- stan_dens(fit, pars="beta[3]") + geom_vline(xintercept = b2) +
  xlab("Q beta")
gridExtra::grid.arrange(plot_b0, plot_b1, plot_b2, plot_phi, plot_sdp, nrow=2)

