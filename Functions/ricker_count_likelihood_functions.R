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
  } else {
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


ricker_count_nb_fit <- function(theta, y, tol = 1e-5, max_iter = 500){
  
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
  } else {
    stop("y must be a vector or an offset dataframe with columns yt and ytm1")
  }
  
  # ---- Objective functions ----
  ## ---- Step 1 ----
  step_1_obj_fun <- function(theta, X, psi){
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
    
    # return the negative log-likelihood
    return(-sum(
      dnbinom(x = X[1:(n - 1), "yt"], mu = exp(eta[2:n]), size = psi, log = T)
    ))
  }
  
  ## ---- Step 2 ----
  step_2_obj_fun <- function(lpsi, X, mu){
    # return the negative log-likelihood
    return(-sum(
      dnbinom(x = X[1:(n - 1), "yt"], mu = mu, size = exp(lpsi), log = T)
    ))
  }
  
  # ---- Main loop ----
  theta_curr <- theta[c("r", "lalpha")]
  lpsi_curr <- theta["lpsi"]
  H <- matrix(nrow = length(theta_curr), ncol = length(theta_curr))
  i <- 0
  diff <- 999
  while(diff > tol & i < max_iter){
    
    # optimize r and lalpha
    fit_i1 <- optim(theta_curr, step_1_obj_fun, X = X, psi = exp(lpsi_curr), hessian = T)
    
    theta_prop <- fit_i1$par
    H <- fit_i1$hessian
    
    # now get mu for the current values
    eta <- vector(mode = "double", length = n)
    eta[1] <- log(X[1, "ytm1"])
    
    # compute means
    for(t in 1:(n - 1)){
      eta[t + 1] <- log(X[t, "ytm1"]) + theta_prop["r"] - exp(theta_prop["lalpha"]) * X[t, "ytm1"]
    }
    
    fit_i2 <- optimize(step_2_obj_fun, interval = c(0, log(50)), X = X, mu = exp(eta))
    
    lpsi_prop <- fit_i2$minimum
    
    diff <- sqrt(sum(
      (c(theta_curr, lpsi_curr) - c(theta_prop, lpsi_prop))^2
    ))
    i <- i + 1
    theta_curr <- theta_prop
    lpsi_curr <- lpsi_prop
    
  }
  
  # ---- compiling return objects ----
  estims <- theta_curr
  estims["lalpha"] <- exp(estims["lalpha"])
  names(estims)[2] <- "alpha"
  
  # CIs and SEs
  V <- solve(H)
  ses <- sqrt(diag(V))
  cis <- mapply(FUN = function(x, s){ x + c(-1, 1) * 2 * s }, x = theta_curr, s = ses)
  # convert lalpha
  cis[, 2] <- exp(cis[, 2])
  cis <- t(cis)
  colnames(cis) <- c("lower_95", "upper_95")
  
  # need to use the delta method for se of alpha
  ses[2] <- sqrt(sqrt(ses[2]) * exp(theta_curr[2])^2)
  names(ses) <- names(estims)
  
  return(
    list(
      estim = estims,
      se = ses,
      lower = as.double(cis[, 1]),
      upper = as.double(cis[, 2]),
      convergence = ifelse(i < max_iter, 0, 1),
      psi = exp(lpsi_curr)
    )
  )
  
}


