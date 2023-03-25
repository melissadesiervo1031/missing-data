# Load packages
library(tidyverse)
library(rstan)
library(forecast)
library(brms)
library(lubridate)



# simulate dataset ####
N <- 500
K <- 2
X <- matrix(c(rep(1, N), rnorm(2*N, 0, 1)), ncol = 3)
beta <- rnorm(K + 1, 0, 2)
sdo <- 1
sdp <- 1

phi <- runif(1, 0, 1)

mu0 <- rnorm(1, X[1,] %*% beta, sdo) 
mu <- rep(mu0, N);

for(i in 2:N){
    mu[i] = X[i, ] %*% beta + mu[i-1] * phi + rnorm(1, 0, sdo)
}

y = rnorm(N, mu, sdp)

# plot(X[,3], mu)
# plot(X[,2], mu)

ARmod <- stan_model("GPP sim and real/Stan_code/AR1_obserror.stan")
ARmod_ss <- stan_model("GPP sim and real/Stan_code/AR1_ss.stan")

datlist <- list(
  N = N, K = K, 
  X = X, y = y)

AR_fit <- rstan::sampling(
  ARmod, datlist,
  chains = 4, cores = 4,
  iter = 4000
  )

ARss_fit <- rstan::sampling(
  ARmod_ss, datlist,
  chains = 4, cores = 4,
  iter = 4000
  )

# compare results of stripped down AR models to the GPP model
AR2 <- stan_model("GPP sim and real/Stan_code/AR1_light_Q_centered.stan")

datlist2 <- list(
  N = N,  
  light = X[,2], 
  Q = X[,3],
  P_obs = y,
  sdo = rep(sdo, N),
  miss_vec = rep(1, N))

AR_f2 <- rstan::sampling(
  AR2, datlist2,
  chains = 4, cores = 4,
  iter = 4000
  )

print(AR_fit, pars = c('phi', 'beta', 'sdo'))
print(ARss_fit, pars = c('phi', 'beta', 'sdo', 'sdp'))
print(AR_f2, pars = c('phi', 'beta', 'sdp'))

# All three of these models correctly recover the parameters, 
# The observation error only model puts all of the error into one term equal 
# to the some of the process and observation error variances. The state space
# model cannot resolve the equifinality between the two error terms, but gets
# the combined error correct. 
# The centered model has the intercept shifted, it can be compared to the other 
# intercepts by multiplying by (1 - phi)

parameters <- list(
  phi = phi,
  beta = beta,
  sdo = sdo, sdp = sdp, sd = sqrt(2))


# compare results to Arima and BRMS models:
simdat <- data.frame(y = y, 
                     l = X[,2],
                     q = X[,3])

arima_mod <- forecast::Arima(mu, order = c(1,0,0), xreg = X[,2:3])
brm_mod <- brms::brm(y ~ ar(p = 1) + l + q, data = simdat)

summary(arima_mod)
summary(brm_mod)
stancode(brm_mod)
# the arima and brms models always seem to recover the same parameters as each other,
# and phi tends to be higher than the true value and than the estimate generated from 
# the Stan models above. Sometimes, the beta parameters are similar to those from the 
# Stan models, sometimes they aren't. Not sure what to make of this.