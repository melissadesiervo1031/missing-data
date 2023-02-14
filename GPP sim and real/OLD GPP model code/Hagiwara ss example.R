library(dlm)

# Introduction and analysis examples of a well-known component model in the linear Gaussian state-space model

## Combination of individual models

## Local-level model

### Example: artificial local-level model

```{r Code 9.1, collapse=TRUE}
# <<Generate artificial data that obey local-level model>>
# Preprocessing
set.seed(23)

# Setting of local-level model
W <- 1
V <- 2
m0 <- 10
C0 <- 9
mod <- dlmModPoly(order = 1, dW = W, dV = V, m0 = m0, C0 = C0)
mode<- dlmModARMA(ar=phi,sigma2=sdp, dV=sdo)
# Generate observations using Kalman prediction
t_max <- 200
sim_data <- dlmForecast(mod = mode, nAhead = t_max, sampleNew = 1)
y <- sim_data$newObs[[1]]
x<- sim_data$newStates[[1]]
# Cast the result to ts class
y <- ts(as.vector(y))
x<- ts(as.vector(x))
# Plot results
plot(y, ylab = "y")
points(x)

sink("model10-1.stan")

cat("

data{
  int<lower=1>   t_max;    // Time series length
  vector[t_max]   y;       // Observations
  
  cov_matrix[1]   W;       // Variance of state noise
  cov_matrix[1]   V;       // Variance of observation noise
  real           m0;       // Mean of prior distribution
  cov_matrix[1]  C0;       // Variance of prior distribution
}

parameters{
  real           x0;       // State [0]
  vector[t_max]   x;       // State [1:t_max]
}

model{
  // Likelihood part
  /* Observation equation; see also equation (5.11) */
    for (t in 1:t_max){
      y[t] ~ normal(x[t], sqrt(V[1, 1]));
    }
  
  // Prior part
  /* Prior distribution for state */
    x0   ~ normal(m0, sqrt(C0[1, 1]));
  
  /* State equation; see also equation (5.10) */
    x[1] ~ normal(x0, sqrt(W[1, 1]));
  for (t in 2:t_max){
    x[t] ~ normal(x[t-1], sqrt(W[1, 1]));
  }
}

"
    ,fill=TRUE)
sink()
closeAllConnections()


stan_mod_out <- stan_model(file = "model10-1.stan")

# Smoothing: execution (sampling)
fit_stan <- sampling(object = stan_mod_out,
                     data = list(t_max = t_max, y = y, 
                                 W = mod$W, V = mod$V, 
                                 m0 = mod$m0, C0 = mod$C0),
                     pars = c("x"),
                     seed = 123
)

oldpar <- par(no.readonly = TRUE); options(max.print = 99999)
fit_stan
par(oldpar)
tmp_tp <- traceplot(fit_stan, pars = c(sprintf("x[%d]", 100), "lp__"), alpha = 0.5)
tmp_tp + theme(aspect.ratio = 3/4)

# Extract necessary sampling results
stan_mcmc_out <- rstan::extract(fit_stan, pars = "x")
str(stan_mcmc_out)
# Calculate the mean, 2.5%, and 97.5% values while marginalizing
s_mcmc <- colMeans(stan_mcmc_out$x)
s_mcmc_quant <- apply(stan_mcmc_out$x, 2, FUN = quantile, probs=c(0.025, 0.975))





sink("model10-2.stan")

cat("


data{
  int<lower=1>   t_max;    // Time series length
  vector[t_max]   y;       // Observations
  
  real           m0;       // Mean of prior distribution
  cov_matrix[1]  C0;       // Variance of prior distribution
}

parameters{
  real           x0;       // State [0]
  vector[t_max]   x;       // State [1:t_max]
  
  cov_matrix[1]   W;       // Variance of state noise
  cov_matrix[1]   V;       // Variance of observation noise
}

model{
  // Likelihood part
  /* Observation equation */
    for (t in 1:t_max){
      y[t] ~ normal(x[t], sqrt(V[1, 1]));
    }
  
  // Prior part
  /* Prior distribution for state */
    x0   ~ normal(m0, sqrt(C0[1, 1]));
  
  /* State equation */
    x[1] ~ normal(x0, sqrt(W[1, 1]));
  for (t in 2:t_max){
    x[t] ~ normal(x[t-1], sqrt(W[1, 1]));
  }
  
  /* Prior distribution for W and V: noninformative prior distribution (utilizing the default setting) */
}

"
    ,fill=TRUE)
sink()
closeAllConnections()

stan_mod_out <- stan_model(file = "model10-2.stan")
# Smoothing: execution (sampling)
fit_stan <- sampling(object = stan_mod_out,
                     data = list(t_max = t_max, y = y, 
                                 m0 = mod$m0, C0 = mod$C0),
                     pars = c("W", "V", "x"),
                     seed = 123
)

traceplot(fit_stan, pars=c("W","V"))

plot_sdo <- stan_dens(fit_stan, pars="W") + geom_vline(xintercept =mod$W)
plot_sdp <- stan_dens(fit_stan, pars="V") + geom_vline(xintercept = mod$V)
grid.arrange(plot_sdp,plot_sdo, nrow=2)
