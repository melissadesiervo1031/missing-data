
###############################################################
# This script generates a bunch of Gaussian AR(1) time series
# datasets with arbitrary covariate matrices and AR parameters
###############################################################


##### Document setup #####


#' Generate random matrix of covariates with MVN distribution
#'
#' @param n Integer number of observations
#' @param p Integer number of covariates (columns)
#' @param sigma2 Positive-valued scalar or vector of variances for the covariates
#' @param mu Scalar or vector of means for the columns
#'
#' @return
#' @export
#'
#' @examples
rand_mod_matrix <- function(n, p, sigma2 = 1, mu = 0, intercept = TRUE){
  M <- matrix(
    data = runif(p^2, min = -1), nrow = p, ncol = p
  )
  # create orthogonal matrix
  L <- qr.Q(qr(M))
  
  # create diagonal matrix of sds or scales
  D <- diag(sigma2, nrow = p, ncol = p)
  
  # create covariance matrix
  Sigma <- t(L) %*% D %*% L
  if(length(mu) == 1){
    mu <- rep(mu, p)
  }
  
  X <- mvtnorm::rmvnorm(
    n = n,
    mean = mu,
    sigma = Sigma
  )
  
  if(intercept){
    return(cbind(rep(1, n), X))
  } else{
    return(X)
  }
  
}




##### Sim parameters #####

# global parameters
  set.seed(235)
  nsims <- 1000
  n <- 365
  p <- 2
  sde <- 1
  
# parameters for each simulation
  params <- purrr::map(
    1:nsims,
    ~ list(
      phi = runif(1, 0, 0.8),
      beta = rnorm(p + 1),
      X = rand_mod_matrix(n = n, p = p)
    )
  )

# simulated datasets
  sims <- purrr::map(
    params,
    ~ list(
      y = {as.double(arima.sim(
        n = n,
        model = list(ar = .x$phi),
        sd = sde
      )) + as.double(.x$X %*% .x$beta)},
      sim_params = .x
    )
  )

# save data
  saveRDS(
    sims,
    file = here::here("data/gauss_ar1_0miss_datasets.rds")
  )


