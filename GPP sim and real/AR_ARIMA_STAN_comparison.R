# Load packages
library(tidyverse)
library(rstan)
library(forecast)
library(brms)
library(lubridate)



# simulate dataset ####
N <- 300
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

plot(X[,3], mu)
plot(X[,2], mu)

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

simdat <- data.frame(y = y, 
                     l = X[,2],
                     q = X[,3])

arima_mod <- forecast::Arima(mu, order = c(1,0,0), xreg = X[,2:3])
brm_mod <- brms::brm(y ~ ar(p = 1) + l + q, data = simdat)

parameters <- list(
  phi = phi,
  beta = beta,
  sdo = sdo, sdp = sdp, sd = sqrt(2))
parameters

print(AR_fit, pars = c('phi', 'beta', 'sdo'))
print(ARss_fit, pars = c('phi', 'beta', 'sdo', 'sdp'))

summary(arima_mod)
summary(brm_mod)

