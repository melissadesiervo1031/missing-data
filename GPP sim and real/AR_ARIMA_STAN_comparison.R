# Load packages
library(tidyverse)
library(rstan)
library(forecast)
library(brms)
library(lubridate)


# compile stan models
# simple AR1 model with observation error only:
ARmod <- stan_model("GPP sim and real/Stan_code/AR1_obserror.stan") 
# state space AR1 model with observation and process error:
ARmod_ss <- stan_model("GPP sim and real/Stan_code/AR1_ss.stan")
# Centered state space AR1 model with observation and process error and missing data loop:
AR2 <- stan_model("GPP sim and real/Stan_code/AR1_light_Q_centered.stan")

# DF <- data.frame()
for(i in 1:10){
  print(paste0("i = ", i))
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


df <- data.frame(
  parameter = c('phi', 'beta1', 'beta2', 'beta3', 'sdo', 'sdp', 'sd'),
  value = c(phi, beta, sdo, sdp, sqrt(2))
)
# plot(X[,3], mu)
# plot(X[,2], mu)

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

get_pars <- function(fit, sdp = FALSE){
  dd <- summary(fit)$summary %>%
    data.frame() %>%slice(1:6) %>%
    select(mean)
  csd = sqrt(dd$mean[5]^2 + dd$mean[6]^2)
  if(!sdp){ dd$mean[6] <- NA_real_
            csd = dd$mean[5]}
  means <- dd$mean[c(4,1,2,3,5,6)]
  means <- c(means, csd)
  
  return(means)
}
df$AR_fit <- get_pars(AR_fit)
df$ARss_fit <- get_pars(ARss_fit, TRUE)
df$AR_f2 <- get_pars(AR_f2, TRUE)


# All three of these models correctly recover the parameters, 
# The observation error only model puts all of the error into one term equal 
# to the some of the process and observation error variances. The state space
# model cannot resolve the equifinality between the two error terms, but gets
# the combined error correct. 
# The centered model has the intercept shifted, it can be compared to the other 
# intercepts by multiplying by (1 - phi)


# compare results to Arima and BRMS models:
simdat <- data.frame(y = y, 
                     l = X[,2],
                     q = X[,3])

arima_mod <- forecast::Arima(mu, order = c(1,0,0), xreg = X[,2:3])
brm_mod <- brms::brm(y ~ ar(p = 1) + l + q, data = simdat)

summary(arima_mod)
df$Arima <- c(arima_mod$coef, rep(NA, 2), arima_mod$sigma2)
a <- summary(brm_mod)
df$brms <- c(a$cor_pars[1,1], a$fixed[,1], rep(NA, 2), a$spec_pars[1,1])

DF <- bind_rows(DF, df)
}
DF <- DF %>% 
  rename(AR_obserr = AR_fit, AR_ss = ARss_fit, AR_missdat = AR_f2) 

write_csv(DF, 'data/AR_model_comparison_test.csv')

DF %>% mutate(across(-parameter, ~(. - value))) %>%
  select(-value) %>%
  pivot_longer(-parameter, names_to = 'model', values_to = 'par_offset') %>%
  ggplot(aes(parameter, par_offset, fill = model)) +
  geom_boxplot() +
  # scale_y_log10() +
  ylim(-5,5)+
  theme_minimal() + geom_hline(yintercept = 0)

# the arima and brms models always seem to recover the same parameters as each other,
# and phi tends to be higher than the true value and than the estimate generated from 
# the Stan models above. Sometimes, the beta parameters are similar to those from the 
# Stan models, sometimes they aren't. Not sure what to make of this.