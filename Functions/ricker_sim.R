


#' Simulate population abundance from a stochastic Ricker model
#'
#' @param n Length of desired time series 
#' @param r Log intrinsic growth factor
#' @param alpha Intra-specific density-dependence
#' @param N0 Initial abundance
#' @param err_fam Error family. Either "poisson" or "neg_binom". If "neg_binom" is 
#' specified, then a dispersion parameter \code{psi} must be included.
#' @param psi Dispersion parameter for negative binomial error distribution.
#'
ricker_sim <- function(n, r, alpha, N0, err_fam = "poisson", psi = NULL){
  
  # initialize vector
  N <- vector(mode = "double", length = n)
  N[1] <- N0
  
  # continue the series
  for(t in 2:n){
    mu_t <- N[t - 1] * exp(r - alpha * N[t - 1])
    if(err_fam == "poisson"){
      N[t] <- rpois(1, lambda = mu_t)
    }
    if(err_fam == "neg_binom"){
      N[t] <- rnbinom(1, mu = mu_t, size = psi)
    }
  }
  
  return(N)
  
}