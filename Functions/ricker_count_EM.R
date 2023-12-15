####################################################
# This is a user-defined functions to fit Ricker
# population models to count data, potentially with 
# missing observations, using EM
####################################################

source(here::here("Functions/ricker_count_likelihood_functions.R"))

#' Fitting the Ricker population model to count data
#' 
#' This function uses an Expectation maximization approach to fit a Ricker population model
#' with Poisson or Negative-Binomial demographic stochasticity. This function can fit the
#'
#' @param y Vector of population counts with NAs placed in positions in which there are missing
#' observations (if any).
#' @param init_theta Initial values for the parameter vector.
#' @param fam Family of the error distribution. Can be either \code{"poisson"} or \code{"neg_binom"}
#' @param tol Tolerance for convergence.
#' @param max_iter Miximum number of iterations to run the EM algorithm before stopping.
#'
#' @return List with estimates and conditional standard errors for \code{theta}, estimates of the
#' full data vector with latent values filled in (\code{z_hat}), and a convergence code. 0 means 
#' the algorithm converged before being stopped.
#' 
ricker_EM <- function(y, init_theta, fam = "poisson", tol = 1e-5, max_iter = 50){
  
  n <- length(y)
  p <- length(init_theta)
  
  # remove starting NAs
  if(is.na(y[1])){
    obs <- which(!is.na(y))
    y <- y[min(obs):n]
    n <- length(y)
  }
  
  # define initial parameter vector
  Theta <- matrix(init_theta, ncol = length(init_theta), nrow= 1)
  
  # define matrix of latent vectors
  Z <- matrix(data = y, ncol = n, nrow = 1)
  
  # iterate through
  dif <- 1
  s <- 1
  
  while(dif > tol & s <= max_iter){
    
    #initialize latent vector
    z_s <- y
    
    # determine parameters sets
    theta <- as.double(Theta[s, ])
    if(fam == "poisson"){
      beta <- theta
    }
    if(fam == "neg_binom"){
      beta <- theta[-p]
      phi <- theta[p]
    }
    
    # fill in missing values with their expected values
    for(t in 2:n){
      if(is.na(z_s[t])){
        z_s[t] <- round(ricker_step(beta, z_s[t - 1]))
        # if rounds to zero, round up instead
        if(z_s[t] == 0){z_s[t] <- 1}
      }
    }
    
    # now find best values of theta given the "full" data
    X <- cbind(
      rep(1, n),
      z_s
    )
    
    # adding try catch here to avoid error "function cannot be evaluated at initial parameters"
    fit_s <- tryCatch({
      optim(
        par = theta, fn = ricker_count_neg_ll,
        y = z_s, X = X, fam = fam, hessian = T
      )
    },error=function(cond){
      message(paste("we have had an error in evaluating function at initial parameters"))
      return(NA)
    })
    
    if(is.na(fit_s[1])){
      return(NA)
    }

    
    # store new value of theta and update
    Theta <- rbind(
      Theta,
      fit_s$par
    )
    Z <- rbind(
      Z, z_s
    )
    
    # calculate the mean relative difference
    if(s == 1){
      dif <- dif
    } else{
      dif <- abs(
        ricker_count_neg_ll(Theta[s, ], Z[s, ], X, fam = fam) - 
          ricker_count_neg_ll(Theta[s + 1, ], Z[s + 1, ], X, fam = fam)
      )
    }
    
    s <- s + 1
    
  }
  
  # compute final values for the return list
  theta_star <- as.double(Theta[nrow(Theta), ])
  convergence <- as.numeric(s == max_iter)
  
  return(list(
    theta = theta_star,
    Theta = Theta,
    z_s = z_s,
    convergence = convergence
  ))
  
}



#' Wrapper to fit a Ricker count model to data using the EM algorithm
#' 
#' This function is a wrapper to fit the stochastic Ricker model with
#' either Poisson or Negative Binomial error distribution and missing observations
#' encoded as NAs.
#'
#' @param y Vector of population counts, with NA in the place of missing observations.
#' @param fam Error family. Can be either c("poisson", "neg_binom").
#' @param ... Additional arguments passed to the EM algorithm, such as initial values 
#' (init_theta = c()) or maximum iterations before stopping (max_iter = 50).
#'
#' @return List of intrinsic growth factor and intra-specific competitive effect estimates,
#' standard errors, and 95% confidence limits. The standard errors are not available using the
#' EM algorithm, so \code{se = NA; lower = NA; upper = NA}. 
#'
#' @examples
#' 
#' y <- readRDS("data/missingDatasets/pois_sim_randMiss_A.rds")[[1]]$y[[1]]
#' fit_ricker_EM(y)
#' 
fit_ricker_EM <- function(y, fam = "poisson", ...){
  
  # Check for population extinction
  if(any(y==0,na.rm=T)){
    warning("population extinction caused a divide by zero problem, returning NA")
    return(list(
      NA,
      cause = "population extinction"
    ))
  }
  
  # Check for NaN
  if(any(is.nan(y),na.rm=T)){
    warning("NaN found, recode missing data as NA, returning NA")
    return(list(
      NA,
      reason = "NaN found"
    ))
  }
  
  # Check for Inf
  if(any(is.infinite(y),na.rm=T)){
    warning("infinite population detected, recheck data returning NA")
    return(list(
      NA,
      reason = "population explosion"
    ))
  }
  
  # remove starting NAs
  if(is.na(y[1])){
    warning("Removing starting NAs...")
    start <- min(which(!is.na(y)))
    y <- y[start:length(y)]
  }
  
  init_theta <- c(0.5, -0.01)
  if(fam == "neg_binom"){
    init_theta <- c(init_theta, 10)
  }
  
  args <- list(
    y = y,
    fam = fam,
    init_theta = init_theta,
    tol = 1e-5, 
    max_iter = 50
  )
  
  args2 <- list(...)
  if(length(args2) > 0){
    args[names(args2)] <- args2
  }
  
  fit <- do.call(ricker_EM, args = args)
  
  # adding condition here to avoid error "function cannot be evaluated at initial parameters"
  if(is.na(fit[1])){
    return(list(NA,reason="we have had an error in evaluating function at initial parameters"))
  }
  
  if(fam == "neg_binom"){
    parnames <- c("r", "alpha", "psi")
    estims <- fit$theta * c(1, -1, 1)
    names(estims) <- parnames
  }
  if(fam == "poisson"){
    parnames <- c("r", "alpha")
    estims <- fit$theta * c(1, -1)
    names(estims) <- parnames
  }
  
  return(list(
    estim = estims,
    se = NA,
    lower = NA,
    upper = NA
  ))
  
}
