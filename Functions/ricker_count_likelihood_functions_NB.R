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
#' @param theta Generic parameter vector. For Poisson: length 2 (r, alpha).
#'   For neg_binom: length 3 (r, alpha, psi) where psi is the size/dispersion
#'   parameter on the natural scale (must be positive).
#' @param y The observed count data
#' @param X (\eqn{n} x \eqn{P}) Model matrix. The first column should be all 1's and the second
#'   should be y. Additional covariates can be added.
#' @param fam Error family: "poisson" or "neg_binom".
#'
#' @return Scalar value of the negative log likelihood of \code{theta} given the data
#'
ricker_count_neg_ll <- function(theta, y, X = NULL, fam = "poisson"){

  n <- length(y)
  p <- length(theta)

  # structural parameters are always the first two (r, alpha)
  theta_struct <- theta[1:2]

  # compute linear predictor
  eta <- vector(mode = "double", length = n)
  eta[1] <- log(y[1])

  if(fam == "poisson"){
    for(t in 2:n){
      eta[t] <- log(y[t - 1]) + X[t - 1, 1:2] %*% theta_struct
    }
    return(-sum(
      dpois(x = y[2:n], lambda = exp(eta[2:n]), log = TRUE)
    ))
  }

  if(fam == "neg_binom"){
    # theta[3] is psi (dispersion/size) on natural scale — must be positive
    psi <- theta[3]
    if(psi <= 0) return(Inf)

    for(t in 2:n){
      eta[t] <- log(y[t - 1]) + X[t - 1, 1:2] %*% theta_struct
    }
    return(-sum(
      dnbinom(x = y[2:n], mu = exp(eta[2:n]), size = psi, log = TRUE)
    ))
  }
}

#' Negative log-likelihood for multiple independent Ricker time series
#'
#' Assumes all series share the same parameter vector \code{theta}. Series are
#' treated as conditionally independent given \code{theta}, so their
#' log-likelihoods are summed.
#'
#' @param theta Parameter vector. For Poisson: c(r, alpha). For neg_binom:
#'   c(r, alpha, psi).
#' @param y_list List of population count vectors (NAs already filled in).
#' @param X_list List of model matrices corresponding to each series.
#' @param fam Error family: "poisson" or "neg_binom".
#'
#' @return Scalar total negative log-likelihood summed across all series.
#'
ricker_count_neg_ll_multi <- function(theta, y_list, X_list, fam = "poisson"){
  total_nll <- 0
  for(i in seq_along(y_list)){
    total_nll <- total_nll + ricker_count_neg_ll(theta, y_list[[i]], X_list[[i]], fam)
  }
  return(total_nll)
}
