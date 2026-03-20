###############################################
# These functions pertain to the log-likelihood
# of a Ricker population model with Poisson or
# Neg-binomial error distribution
###############################################

#' One step in the deterministic Ricker population process
#'
#' @param theta Parameter vector with intrinsic growth rate \eqn{r} as the first element
#' and negative density dependence \eqn{alpha} as the second parameter.
#' NOTE that \eqn{alpha} should be negative to induce negative density dependence.
#' @param Nt The population size at the current time.
#'
#' @return Population size at the next time step.
#' 
ricker_step <- function(theta, Nt){
  return(
    Nt * exp(theta[1] + theta[2] * Nt)
  )
}




#' Negative log-likelihood function for Ricker time series with
#' conditionally Poisson- or negative-binomial-distributed observations
#'
#' @param theta Generic parameter vector of length \eqn{P}.
#' @param y The observed count data
#' @param X (\eqn{n} x \eqn{P}) Model matrix. The first column should be all 1's and the second
#' should be y. Additional covariates can be added.
#'
#' @return Scalar value of the negative log likelihood of \code{theta} given the data
#' 
ricker_count_neg_ll <- function(theta, y, X = NULL, fam = "poisson"){
  
  n <- length(y)
  p <- length(theta)
  # compute means
  eta <- vector(mode = "double", length = n)
  eta[1] <- log(y[1])
  
  if(fam == "poisson"){
    for(t in 2:n){
      eta[t] <- log(y[t - 1]) + X[t - 1, ] %*% theta
    }
    
    # return the negative log-likelihood
    return(-sum(
      dpois(x = y[2:n], lambda = exp(eta[2:n]), log = T)
    ))
  }
  if(fam == "neg_binom"){
    b <- theta[1:(p-1)]
    for(t in 2:n){
      eta[t] <- log(y[t - 1]) + X[t - 1, ] %*% b
    }
    
    # return the negative log-likelihood
    return(-sum(
      dnbinom(x = y[2:n], mu = exp(eta[2:n]), size = theta[p], log = T)
    ))
  }
  
}


zt_poisson_rng <- function(n, lambda){

  p0 <- dpois(0, lambda = lambda)
  u <- runif(n, min = p0, max = 1)
  return(
    qpois(u, lambda = lambda)
  )

}


zt_neg_binom_rng <- function(n, mu, size){
    p0 <- dnbinom(0, size = size, mu = mu)                                      
    u  <- runif(n, min = p0, max = 1)                                           
    return(qnbinom(u, size = size, mu = mu))                                    
}


#' Negative log-likelihood for multiple independent Ricker time series
#' 
#' @param theta Parameter vector
#' @param y_list List of population count vectors (NAs already filled in)
#' @param X_list List of model matrices corresponding to each series
#' @param fam Error family
#' 
ricker_count_neg_ll_multi <- function(theta, y_list, X_list, fam = "poisson"){
  total_nll <- 0
  for(i in seq_along(y_list)){
    total_nll <- total_nll + ricker_count_neg_ll(theta, y_list[[i]], X_list[[i]], fam)
  }
  return(total_nll)
}


