
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

# Simulate data
# simulate arima process

mod_dat <- data.frame()
for(i in 1:100){
if(i %% 10 == 0) print(paste0('model ', i))
n <- 365
p <- 2
phi <- runif(1, 0, 0.8) ###phi gets weird too close to 1 ###
beta <- rnorm(p + 1)
X <- matrix(c(rep(1, n), rnorm(n = n * p)), nrow = n, ncol = p+1)
sde <- 1
mu <- as.double(X %*% beta)

y <- arima.sim(
   n = n,
   model = list(ar = phi),
   sd = sde) + mu

#### Parameter recovery using arima and brms
# phi; beta; sde
fit <- arima(y, order = c(1,0,0), xreg = X[,-1])

dat <- data.frame(y = y,
                  x1 = X[,2],
                  x2 = X[,3])

bform <- brms::bf(y ~ x1 + x2 + ar(p = 1))
bmod <- brms::brm(bform, data = dat)

# brms::posterior_interval(bmod)

#### Comparison to data augmentation in STAN ###

datlist <- list(N = n,
                y = y, 
                light = X[,2],
                discharge = X[,3])

ar_fit <- stan('GPP sim and real/Stan_code/AR1.stan',
                data = datlist,
                chains = 4,
                iter = 2000)
ar_fit2 <- stan('GPP sim and real/Stan_code/AR1_bob.stan',
                data = datlist,
                chains = 4,
                iter = 2000)

# phi; beta; sde
# fit
# bmod
# print(ar_fit, pars = c('phi', 'beta', 'sigma'))

dd <- data.frame(pars = c('phi', 'beta1', 'beta2', 'beta3', 'sigma'), 
                 value = c(phi, beta, sde),
                 arima = c(fit$coef, fit$sigma2), 
                 brms = c(summary(bmod)$cor_pars[,1], summary(bmod)$fixed[,1], 
                          summary(bmod)$spec_pars[,1]),
                 stan = c(summary(ar_fit)$summary[c(4, 1,2,3, 5), 1]),
                 bob = c(summary(ar_fit2)$summary[c(4, 1,2,3,5), 1]),
                 model = paste0('mod_', i))

row.names(dd) <- NULL

mod_dat <- bind_rows(mod_dat, dd)
}


md <- mod_dat %>%
    # group_by(model) %>%
    # mutate(stan = case_when(pars == 'beta1' ~ stan/(1 - stan[pars == 'phi']),
    #                         TRUE ~ stan)) %>%
    mutate(arima_diff = arima - value,
           brm_diff = brms - value,
           stan_diff = stan - value, 
           bob_diff = bob - value) 
  

md %>%
    pivot_longer(cols = ends_with('diff'), 
                 values_to = 'parameter_difference',
                 names_to = 'fit') %>%
    ggplot(aes(pars, parameter_difference, fill = fit)) +
    geom_boxplot()+
    theme_classic()

md %>% 
  ggplot(aes(value, bob)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(.~pars) +
  theme_classic()

md %>% 
  ggplot(aes(value, arima)) +
  geom_point() +
  geom_point(aes(y = stan), col = 'brown') +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(.~pars) +
  theme_classic()

md %>% 
    group_by(model) %>%
    mutate(phi_value = value[pars == 'phi'])%>%
    pivot_longer(cols = ends_with('diff'), 
                 values_to = 'parameter_difference',
                 names_to = 'fit') %>% 
    pivot_wider(id_cols = c('model', 'fit', 'phi_value'), 
                values_from = 'parameter_difference',
                names_from = 'pars') %>%
    ggplot(aes(phi_value, phi, col = fit))+
    geom_point() + geom_abline(intercept = 0, slope = 0)+
    ylab('Error in phi estimate') +
    xlab('phi value') + theme_classic()

ggpubr::ggarrange(a,b, ncol = 1, common.legend = TRUE)

