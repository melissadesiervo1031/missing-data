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
    b <- theta[names(theta) != "psi"]
    for(t in 2:n){
      eta[t] <- log(y[t - 1]) + X[t - 1, ] %*% b
    }
    
    # return the negative log-likelihood
    return(-sum(
      dnbinom(x = y[2:n], mu = exp(eta[2:n]), size = theta["psi"], log = T)
    ))
  }
  
}


zt_poisson_rng <- function(n, lambda){
  
  u <- runif(n)
  
  t <- -log(1 - u * (1 - exp(-lambda)))
  
  lambda_prime <- lambda - t
  
  return(rpois(n, lambda_prime) + 1)
  
}


zt_neg_binom_rng <- function(n, mu, size, max_iter = 1000){
  y <- rnbinom(n, size = size, mu = mu)
  i <- 1
  while(sum(y == 0) > 0 & i < max_iter){
    n_new <- sum(y == 0)
    y_new <- rnbinom(n_new, size = size, mu = mu)
    y[y == 0] <- y_new
    i <- i + 1
  }
  if(i == max_iter){
    stop("Timed out. Increase the max iterations or try an alternative approach.")
  }
  return(y)
}

