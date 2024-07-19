# Simulate GPP data from AR1_light_Q_centered model
# Author: AM Carter


# load packages -----------------------------------------------------------
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())


# load data ---------------------------------------------------------------
# load data from NWIS
dat <- read_csv('data/NWIS_MissingTS_subset.csv')
# load metadata from NWIS
mdat <- read_csv('data/NWIS_MissingTSinfo_subset.csv')


# select a single site
id <- mdat$site_name[1]
dd <- filter(dat, site_name == id)

# calculate observation error based on 95% CIs

sigma_obs <- (dd$GPP.upper - dd$GPP.lower)/3.92
dd <- dd %>%
  select(date, GPP, light, Q) %>%
  mutate(across(-date, ~zoo::na.approx(.)))%>%
  mutate(light.rel = light/max(light),
         Q = scale(Q)[,1])

#Set Parameter values: ####
n <- nrow(dd)
phi <- 0.6
beta <- c(1, 3, -0.5)
sdp <- 1

# randomly assign a chunk of missing data:
miss_vec <- sample(c(rep(1, 320), rep(0, n - 320)), n, replace = FALSE) 
mu <- vector()
mu[1] <- rnorm(1, dd$GPP[1], sdp)

for(i in 2:n){
    mu[i] = rnorm(1, beta[1]*(1-phi) + mu[i-1]*phi + dd$light.rel[i]*beta[2] + dd$Q[i]*beta[3], sdp);
}

P_obs <- rnorm(n, mu, sigma_obs)

model_lq <- stan_model("GPP sim and real/Stan_code/AR1_light_Q_centered.stan")

#Create data object
data <- list(N = n, P_obs = P_obs,
             sdo = sigma_obs, light = dd$light.rel, Q = dd$Q,
             miss_vec = rep(1, nrow(pr)))

#Run Stan
fit_lq <- rstan::sampling(object=model_lq, data = data,  
                          iter = 4000, chains = 4)

# fit with missing data
data <- list(N = n, P_obs = P_obs,
             sdo = sigma_obs, light = dd$light.rel, Q = dd$Q,
             miss_vec = miss_vec)
fit_lq_miss <- rstan::sampling(object=model_lq, data = data,  
                          iter = 4000, chains = 4)

# examine model outputs
traceplot(fit_lq, pars=c("phi", "sdp", "beta"))
pairs(fit_lq, pars=c("phi", "sdp","beta","lp__"))
plot(fit_lq, pars=c("phi", "sdp", "beta"))
print(fit_lq, pars=c("phi", "sdp", "beta"))

traceplot(fit_lq_miss, pars=c("phi", "sdp", "beta"))
pairs(fit_lq_miss, pars=c("phi", "sdp","beta","lp__"))
plot(fit_lq_miss, pars=c("phi", "sdp", "beta"))
print(fit_lq_miss, pars=c("phi", "sdp", "beta"))


# look at posterior predictions:
get_pp <- function(fit){
  fit_extract <- rstan::extract(fit)
  pp <- t(apply(fit_extract$GPP_rep, 2, 
                function(x) quantile(x, probs = c(0.025, 0.5, 0.975), names = F))) %>%
    data.frame() 
  names(pp) <- c('GPP_rep.lower', 'GPP_rep', 'GPP_rep.upper')
  
  return(pp)
  
}

dd$GPP <- P_obs
post <- bind_cols(dd, get_pp(fit_lq))
plot_gpp(post, pp_fit = TRUE)

dd$GPP[miss_vec == 0] <- NA_real_
post_miss <- bind_cols(dd, get_pp(fit_lq_miss))
plot_gpp(post_miss, pp_fit = TRUE)
