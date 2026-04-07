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
    Nt * exp(theta[1] - theta[2] * Nt)
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



#' Constrained negative log-likelihood for a Ricker count time series
#'
#' Computes the negative log-likelihood of a Ricker population model with
#' Poisson or Negative Binomial observations, where the density-dependence
#' parameter \eqn{\alpha} is log-parameterized to enforce positivity. The
#' log-linear mean at time \eqn{t} is:
#' \deqn{\log \mu_t = \log N_{t-1} + r - \alpha N_{t-1}}
#' where \eqn{\alpha = \exp(\texttt{lalpha}) > 0}, ensuring negative density
#' dependence. This is the likelihood used by fitting routines that optimize
#' over a constrained parameter space.
#'
#' Unlike \code{\link{ricker_count_neg_ll}}, this function accepts \code{y}
#' in several formats and constructs the lag pairs internally, so no separate
#' model matrix needs to be supplied.
#'
#' @param theta Named numeric vector of parameters. Must contain:
#'   \describe{
#'     \item{\code{r}}{Intrinsic growth rate (unconstrained, real-valued).}
#'     \item{\code{lalpha}}{Log of the intraspecific competition coefficient.
#'       Exponentiated internally so that \eqn{\alpha = e^{\texttt{lalpha}} > 0}.}
#'     \item{\code{lpsi}}{Log of the Negative Binomial dispersion parameter
#'       (\eqn{\psi = e^{\texttt{lpsi}} > 0}). Required when \code{fam = "neg_binom"};
#'       ignored otherwise.}
#'   }
#' @param y Observed count data. Accepted formats:
#'   \describe{
#'     \item{Numeric vector}{Raw time series. Lag pairs \code{(yt, ytm1)} are
#'       constructed internally.}
#'     \item{Data frame or matrix}{Must have named columns \code{yt} (count at
#'       time \eqn{t}) and \code{ytm1} (count at time \eqn{t-1}), i.e., the
#'       offset/lag-pair format used elsewhere in this package. Rows may come
#'       from multiple concatenated series when using \code{off_patch} mode.}
#'   }
#'   \code{NA}s should be filled in before calling this function.
#' @param X Ignored. Retained for interface compatibility with
#'   \code{\link{ricker_count_neg_ll}}. Lag pairs are always derived from \code{y}.
#' @param fam Error distribution family. One of \code{"poisson"} (default) or
#'   \code{"neg_binom"}.
#'
#' @return Scalar negative log-likelihood evaluated at \code{theta}.
#'
#' @seealso \code{\link{ricker_count_neg_ll}} for the unconstrained version,
#'   \code{\link{ricker_step}} for the deterministic Ricker step,
#'   \code{\link{fit_ricker_cc}} and \code{\link{fit_ricker_EM}} which use
#'   this function internally.
#'
ricker_count_neg_ll_cnstr <- function(theta, y, X = NULL, fam = "poisson"){
  
  if(is.vector(y)){
    X <- data.frame(
      yt = y[2:length(y)],
      ytm1 = y[1:(length(y) - 1)]
    )
  } else if(is.data.frame(y)){
    if(!all(c("yt", "ytm1") %in% names(y))){
      stop("Expecting an offset dataframe with columns yt and ytm1")
    }
    X <- y[, c("yt", "ytm1")]
  } else if(is.matrix(y)){
    if(!all(c("yt", "ytm1") %in% colnames(y))){
      stop("Expecting an offset dataframe with columns yt and ytm1")
    }
    X <- y[, c("yt", "ytm1")]
  } else{
    stop("y must be a vector or an offset dataframe with columns yt and ytm1")
  }
  
  n <- nrow(X) + 1
  p <- length(theta)
  
  # constrained alpha
  alpha <- exp(theta["lalpha"])
  
  # compute means
  eta <- vector(mode = "double", length = n)
  eta[1] <- log(X[1, "ytm1"])

  # compute means
  for(t in 1:(n - 1)){
    eta[t + 1] <- log(X[t, "ytm1"]) + theta["r"] - alpha * X[t, "ytm1"]
  }
  
  if(fam == "poisson"){
    
    # return the negative log-likelihood
    return(-sum(
      dpois(x = X[1:(n - 1), "yt"], lambda = exp(eta[2:n]), log = T)
    ))
  }
  if(fam == "neg_binom"){
    psi <- exp(theta["lpsi"])

    # return the negative log-likelihood
    return(-sum(
      dnbinom(x = X[1:(n - 1), "yt"], mu = exp(eta[2:n]), size = psi, log = T)
    ))
  }
  
}


#' Draw samples from a zero-truncated Poisson distribution
#'
#' Generates random draws from a Poisson distribution conditioned on the outcome
#' being strictly positive (\eqn{X \geq 1}). Uses the inverse-CDF method: a
#' uniform variate is drawn on \eqn{[p_0, 1]}, where \eqn{p_0 = P(X = 0)} under
#' the untruncated Poisson, and then mapped through the Poisson quantile function.
#' This guarantees all returned values are \eqn{\geq 1} without rejection sampling.
#' Used during data augmentation to impute missing counts while avoiding zeros
#' that would cause \eqn{\log(0)} in the likelihood.
#'
#' @param n Number of draws to return.
#' @param lambda Rate parameter of the underlying (untruncated) Poisson distribution.
#'   Must be positive.
#'
#' @return Integer vector of length \code{n} with all elements \eqn{\geq 1}.
#'
#' @seealso \code{\link{zt_neg_binom_rng}} for the Negative Binomial analogue.
#'
zt_poisson_rng <- function(n, lambda){

  p0 <- dpois(0, lambda = lambda)
  u <- runif(n, min = p0, max = 1)
  return(
    qpois(u, lambda = lambda)
  )

}


#' Draw samples from a zero-truncated Negative Binomial distribution
#'
#' Generates random draws from a Negative Binomial distribution conditioned on
#' the outcome being strictly positive (\eqn{X \geq 1}). Uses the same
#' inverse-CDF method as \code{\link{zt_poisson_rng}}: a uniform variate is
#' drawn on \eqn{[p_0, 1]}, where \eqn{p_0 = P(X = 0)} under the untruncated
#' Negative Binomial, and then mapped through the NB quantile function.
#' Used during data augmentation to impute missing counts under Negative
#' Binomial error while avoiding zeros in the likelihood.
#'
#' @param n Number of draws to return.
#' @param mu Mean of the underlying (untruncated) Negative Binomial distribution.
#'   Must be positive. Corresponds to the \code{mu} parameterization in
#'   \code{\link[stats]{dnbinom}}.
#' @param size Dispersion (size) parameter of the Negative Binomial distribution
#'   (\eqn{\psi} elsewhere in this package). Must be positive. Larger values
#'   approach the Poisson limit.
#'
#' @return Integer vector of length \code{n} with all elements \eqn{\geq 1}.
#'
#' @seealso \code{\link{zt_poisson_rng}} for the Poisson analogue,
#'   \code{\link[stats]{dnbinom}} for the parameterization used.
#'
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



