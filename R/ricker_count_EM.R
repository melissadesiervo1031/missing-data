####################################################
# This is a user-defined functions to fit Ricker
# population models to count data, potentially with 
# missing observations, using EM
####################################################



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
      beta <- theta[p]
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
  V <- solve(fit_s$hessian)
  convergence <- as.numeric(s == max_iter)
  
  return(list(
    theta = theta_star,
    Theta = Theta,
    V_c = V,
    z_s = z_s,
    convergence = convergence
  ))
  
}


