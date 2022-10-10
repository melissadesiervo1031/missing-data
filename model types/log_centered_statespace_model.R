# setwd('C:/Users/alice.carter/git/missing-data/')

library(rstan)
library(lubridate)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#Simulate data

N <- 365 #length of data
z <- numeric(N + 1)
tt <- 1:N
time <- seq(as.POSIXct("2020/1/1"),as.POSIXct("2020/12/30"), by="day") #for light


##light
# From Yard et al. (1995) Ecological Modelling.  Remember your trig?
# calculate light as umol photon m-2 s-1.
# Arguments are:
# time = a date and time input (posixct object)
# lat = latitude of field site
# longobs = longitude of field site
# longstd = standard longitude of the field site (NOTE: watch daylight savings time!!!). For PST, longstd is be 120 degrees. But during PDT it is 105 degrees. MST is 105 deg. MDT is 90.

# convert degrees to radians
radi <- function(degrees){(degrees*pi/180)}

# function to estimate light
lightest<- function (time, lat, longobs, longstd) {
  jday<-yday(time)
  E<- 9.87*sin(radi((720*(jday-81))/365)) - 7.53*cos(radi((360*(jday-81))/365)) - 1.5*sin(radi((360*(jday-81))/365))
  LST<-as.numeric(time-trunc(time))
  ST<-LST+(3.989/1440)*(longstd-longobs)+E/1440
  solardel<- 23.439*sin(radi(360*((283+jday)/365)))
  hourangle<-(0.5-ST)*360
  theta<- acos( sin(radi(solardel)) * sin(radi(lat)) + cos(radi(solardel)) * cos(radi(lat)) * cos(radi(hourangle)) )
  suncos<-ifelse(cos(theta)<0, 0, cos(theta))
  GI<- suncos*2326
  GI

}

#Estimate 'real' light data
light<-lightest(time, 47.8762, -114.03, 105) #Flathead Bio Station just for fun

#Condense light to be relative to max light
light.rel<-light/max(light)

##GPP based on time-series model with known parameters
# GPP_t = b0 + phi*GPP_t-1 + b1*L + proc_error
# GPP_obs = GPP + obs_error
b0 <- 1      # intercept
b1 <- 0.6   # light coefficient
phi <- 0.4   # ar1 coefficient
sigma_proc <- 0.1 # ~10% error
sigma_obs <- 0.01

# Expected value based on parameters:  ~10.4
EV <- exp((b0*(1-phi) + b1 * light.rel[1])/(1-phi))
logP <- numeric(N)  # GPP state variable
logP[1] <- log(EV)

# Observed GPP, on a log scale
logP_obs = rnorm(1, logP[1], sigma_obs)
for (i in 2:365){
  logP[i] <- b0*(1-phi) + b1 * light.rel[i] + phi * logP[i-1] + rnorm(1, 0, sigma_proc)
  logP_obs[i] = rnorm(1, logP[i], sigma_obs)
}

plot(exp(logP), type = 'l', ylab = 'GPP')
points(exp(logP_obs))
legend('topright', c('state', 'observations'), pch = c(NA, 1), lty = c(1, 0),
       bty = 'n', ncol = 2)

##Bind simulated data and time
y_full <- as.data.frame(cbind(z = logP, y = logP_obs, tt))
y_full$tt <- as.factor(y_full$tt)

#Stan model
sink("ss_light_log_cent.stan")

cat("

/*----------------------- Data --------------------------*/
  data {
    int N; // Length of state and observation time series
    vector[N] logP_obs; // log observations
    vector[N] light;
    real logP_0; //initial underlying state
  }

/*----------------------- Parameters --------------------------*/

  parameters {
    vector <lower = 0> [N] logP; //log latent state variable
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
    sdo ~ normal(0.01,0.001);

    // Distribution for the zero and first state
    logP[1] ~ normal(logP_0,sdp);

    // Distributions for states
    for(i in 2:N){
       logP[i] ~ normal(b0*(1-phi) + logP[i-1]*phi + light[i]*b1, sdp);// process model with error
        }

    // Observation error model
    logP_obs ~ normal(logP, sdo);

      }" ,fill=TRUE)

sink()
closeAllConnections()

#Prep model
model<-"ss_light_log_cent.stan"
model<-stan_model(model)

#Create data object
data <- list(N = length(logP_obs), logP_obs = logP_obs,
             logP_0 = logP_obs[1], light=light.rel)

#Run Stan
fit <- rstan::sampling(object=model, data = data,  iter = 2000, chains = 4)

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
