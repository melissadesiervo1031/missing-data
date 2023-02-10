####################################################
# These are all user-defined functions that pertain
# to simulating and fitting Ricker population models
# to data, potentially with missing observations
####################################################



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
#' condtionally Poisson-distributed observations
#'
#' @param theta Generic parameter vector of length \eqn{P}.
#' @param y The observed count data
#' @param X (\eqn{n} x \eqn{P}) Model matrix. The first column should be all 1's and the second
#' should be y. Additional covariates can be added.
#'
#' @return Scalar value of the negative log likelihood of \code{theta} given the data
#' 
ricker_count_neg_ll <- function(theta, y, X, fam = "poisson"){
  
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
    for(t in 2:n){
      eta[t] <- log(y[t - 1]) + X[t - 1, ] %*% theta[-p]
    }
    
    # return the negative log-likelihood
    return(-sum(
      dnbinom(x = y[2:n], mu = exp(eta[2:n]), size = theta[p], log = T)
    ))
  }
  
}




suppEM <- function(y, Theta, fam = "poisson"){
  
  
  
}





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
ricker_fit_EM <- function(y, init_theta, fam = "poisson", tol = 1e-5, max_iter = 50){
  
  n <- length(y)
  p <- length(init_theta)
  if(init_theta[1] > 0.7 | init_theta[2] > 0){
    warning("For biologically reasonable starting values for r and alpha, 
            set -0.7 < init_theta[1] < 0.7 and -0.1 < init_theta[2] < 0. 
            The algorithm may not converge with these starting values.")
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
      beta <- theta
      phi <- theta[p]
    }
    
    # fill in missing values with their expected values
    for(t in 2:n){
      if(is.na(z_s[t])){
        z_s[t] <- round(ricker_step(beta, z_s[t - 1]))
      }
    }
    
    # now find best values of theta given the "full" data
    X <- cbind(
      rep(1, n),
      z_s
    )
    fit_s <- optim(
      par = theta, fn = ricker_count_neg_ll,
      y = z_s, X = X, fam = fam, hessian = T
    )
    
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
  
  # if the algorithm has converged, compute standard errors
  if(convergence == 0){
    if(sum(is.na(y)) == 0){
      theta_ses <- sqrt(diag(solve(fit_s$hessian)))
    }
    # adjust for missing information if any
    if(sum(is.na(y)) > 0){
      DM <- 
    }
  }
  
  return(list(
    theta = theta_star,
    wald_se = theta_ses,
    z_hat = z_s,
    convergence = convergence
  ))
  
}

