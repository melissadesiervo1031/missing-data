
###############################################################
# This script generates a bunch of integer-valued time series
# datasets using a Ricker population model
###############################################################


ricker_sim <- function(n, r, alpha, N0){
  
  # initialize vector
  N <- vector(mode = "double", length = n)
  N[1] <- N0
  
  # continue the series
  for(t in 2:n){
    mu_t <- N[t - 1] * exp(r - alpha * N[t - 1])
    N[t] <- rpois(1, lambda = mu_t)
  }
  
  return(N)
  
}



##### Parameters #####

# global parameters
  set.seed(5681)
  nsims <- 1000
  n <- 60
  
# per-simulation parameters
  params <- purrr::map(
    1:nsims,
    ~ list(
      r = runif(1, max = log(2)),
      alpha = runif(1, min = 0.001, 0.01),
      N0 = rpois(1, lambda = 20)
    )
  )
  
##### Simulate the data #####
  sims <- purrr::map(
    params,
    ~ list(
      y = ricker_sim(
        n = n,
        r = .x$r,
        alpha = .x$alpha,
        N0 = .x$N0
      ),
      sim_params = .x
    )
  )
  
  
##### Save data #####
  saveRDS(
    sims,
    file = here::here("data/ricker_0miss_datasets.rds")
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
