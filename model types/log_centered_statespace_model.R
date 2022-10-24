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
b0 <- 1      # intercept
b1 <- 0.6   # light coefficient
phi <- 0.4   # ar1 coefficient
sigma_proc <- 0.1 # ~10% error
sigma_obs <- 0.1

# Expected value based on parameters:  ~10.4
EV <- exp((b0*(1-phi) + b1 * light.rel[1])/(1-phi))
GPP <- numeric(N)  # GPP state variable
GPP[1] <- EV

# Observed GPP: add observation error on a linear, not a log scale
P_obs = rnorm(1, GPP[1], sigma_obs)
for (i in 2:N){
  GPP[i] <- exp(b0*(1-phi) + b1 * light.rel[i] + phi * log(GPP[i-1]) + 
                  rnorm(1, 0, sigma_proc))
  P_obs[i] = rnorm(1, GPP[i], sigma_obs)
}

plot(GPP, type = 'l', ylab = 'GPP')
points(P_obs)
legend('topright', c('state', 'observations'), pch = c(NA, 1), lty = c(1, 0),
       bty = 'n', ncol = 2)

##Bind simulated data and time
y_full <- as.data.frame(cbind(GPP, P_obs, doy = as.factor(tt)))

#Stan model
sink("ss_light_log_cent.stan")

cat("

/*----------------------- Data --------------------------*/
  data {
    int N; // Length of state and observation time series
    vector[N] P_obs; // observed GPP
    vector[N] light;
    real logP_0; // log initial value for GPP
  }

/*----------------------- Parameters --------------------------*/

  parameters {
    vector <lower = 0> [N] logGPP; // log latent variable - true GPP
    real<lower=0> sdp; // Standard deviation of the process equation
    real<lower=0> sdo; // Standard deviation of the observation equation
    real b0; // intercept
    real b1; // light coefficient
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
  }

/*----------------------- Model --------------------------*/

  model {
    // Prior distributions
    phi ~ beta(1,1);
    b0 ~ normal(0, 5);
    b1 ~ normal(0, 5);
    sdp ~ normal(0,1);
    sdo ~ normal(0,0.5);

    // Distribution for the zero and first state
    logGPP[1] ~ normal(logP_0, sdp);

    // Distributions for states
    for(i in 2:N){
       logGPP[i] ~ normal(b0*(1-phi) + logGPP[i-1]*phi + light[i]*b1, sdp); // process model with error
        }

    // Observation error model
    P_obs ~ normal(exp(logGPP), sdo);

      }
    
    " ,fill=TRUE)

sink()
closeAllConnections()

#Prep model
model <- "ss_light_log_cent.stan"
model <- stan_model(model)

#Create data object
data <- list(N = length(P_obs), P_obs = P_obs,
             logP_0 = log(P_obs[1]), light=light.rel)

#Run Stan
fit <- rstan::sampling(object=model, data = data,  
                       iter = 4000, chains = 4)

##Print and extract
fit_extract <- rstan::extract(fit)
print(fit, pars=c( "sdo","phi", "sdp", "b0","b1" ))
pairs(fit, pars=c("sdo","phi", "sdp","b0","b1","lp__"))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit)

##Traceplots
traceplot(fit, pars=c("sdo","phi", "sdp", "b1", "b0"))

##Plot density plots
plot_sdo <- stan_dens(fit, pars="sdo") + geom_vline(xintercept = sigma_obs)
plot_sdp <- stan_dens(fit, pars="sdp") + geom_vline(xintercept = sigma_proc)
plot_phi <- stan_dens(fit, pars="phi") + geom_vline(xintercept = phi)
plot_b0 <- stan_dens(fit, pars="b0") + geom_vline(xintercept = b0) +
  xlab("Intercept")
plot_b1 <- stan_dens(fit, pars="b1") + geom_vline(xintercept = b1) +
  xlab("Light beta")
gridExtra::grid.arrange(plot_b0, plot_b1, plot_phi, plot_sdp, plot_sdo, nrow=2)

# I have found that recovering observation error is really difficult in a log
# model, no matter what value is chosen. I think that while incorporating this
# error structure does provide a more realistic description of the data, it
# probably isn't the move for the missing data paper.
