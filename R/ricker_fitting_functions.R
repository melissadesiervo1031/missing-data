####################################################
# These are all user-defined functions that pertain
# to fitting Ricker population models to count data, 
# potentially with missing observations
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







#' Supplemental EM aglorithm for estimating the covariance matrix of parameter estimates
#' 
#' This function implements the SEM algorithm of Meng and Rubin (1991) in order to esimate
#' the covariance of the parameter estimates. This allows for standard error estimates that
#' account for the amount of missing information.
#'
#' @param y Vector of integer-valued responses with NA in place of missing data
#' @param Theta (S x p) matrix of parameter values from the EM output, where S is the
#' number of steps before convergence for the EM algorithm and p is the number of parameters.
#' @param theta_hat MLE of the parameter vector.
#' @param V_c Complete covariance matrix obtained by inverting the Hessian of the log-likelihood
#' evaluated at the MLE, assuming no missing data (i.e., treating the expected values of the missing
#' data as observed data).
#' @param fam Error distribution. Can be either \code{"poisson"} or \code{"neg_binom"}.
#' @param tol Numerical error tolerance for convergence.
#' @param max_iter Maximum number of iterations over which to run the algorithm.
#'
#' @return Estimated covariance matrix for the parameter estimates.
#'
suppEM <- function(y, Theta, V_c, fam = "poisson", tol = 1e-5, max_iter = 50){
  
  # useful variables
  n <- length(y)
  theta_hat <- Theta[nrow(Theta), ]
  p <- length(theta_hat)
  t <- max(nrow(Theta) - 2, 1)
  
  # initialize matrix of interest
  R <- matrix(data = NA, nrow = p, ncol = p)
  
  # convergence codes
  convergence <- vector(mode = "double", length = p)
  
  # loop through rows to find r_ij-star
  for(i in 1:p){
    
    # initialize inputs for SEM algorithm
    theta_is <- theta_hat
    theta_is[i] <- Theta[t, i]
    theta_isp1 <- ricker_EM(y, init_theta = theta_is, fam = fam, max_iter = 1)$theta
    r_is <- (theta_isp1 - theta_hat) / (theta_is[i] - theta_hat[i])
    
    # loop to find stable r_ij
    dif <- 1; s <- 1
    while(dif > tol & s <= max_iter){
      
      theta_is <- theta_isp1
      theta_isp1 <- ricker_EM(y, init_theta = theta_is, fam = fam, max_iter = 1)$theta
      r_isp1 <- (theta_isp1 - theta_hat) / (theta_is[i] - theta_hat[i])
      dif <- sum(abs(r_is - r_isp1))
      
      # recycle objects
      r_is <- r_isp1
      s <- s + 1
      
    }
    
    R[i, ] <- r_is
    
    # did we converge?
    convergence[i] <- as.numeric(s == max_iter)
    
  }
  
  return(
    V_c %*% (diag(nrow = p) - R)
  )
  
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
ricker_EM <- function(y, init_theta, fam = "poisson", tol = 1e-5, max_iter = 50){
  
  n <- length(y)
  p <- length(init_theta)
  
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









fit_ricker <- function(y, fam, X = NULL, init_theta = NULL, tol = 1e-5, max_iter = 50){
  
  # initial values for parameter vector if not supplied by the user
  if(is.null(init_theta)){
    if(is.null(X)){px <- 0} else{px <- ncol(X)}
    if(fam == "neg_binom"){px <- px + 1}
    p <- 2 + px
    
    # now initialize some parameter values
    init_theta[1] <- runif(1, max = 0.7)
    init_theta[2] <- runif(1, min = -0.1, max = 0)
    if(fam == "neg_binom"){
      init_theta[3:(p - 1)] <- rnorm(p - 3)
      init_theta[p] <- rexp(1)
    }
    if(fam == "poisson"){
      init_theta[3:p] <- rnorm(p - 2)
    }
    
  }
  
  # warning message about biologically reasonable starting values
  if(init_theta[1] > 0.7 | init_theta[2] > 0){
    warning("For biologically reasonable starting values for r and alpha, 
            set 0 < init_theta[1] < 0.7 and -0.1 < init_theta[2] < 0. 
            The algorithm may not converge with these starting values.")
  }
  
  # useful variables
  n <- length(y)
  
  # construct model matrix
  if(!is.null(X)){
    X <- cbind(
      rep(1, n), y,
      X
    )
  } else{
    X <- cbind(
      rep(1, n), y
    )
  }
  
  # construct names for parameter
  pnames <- c("r", "alpha")
  if(fam == "neg_binom"){
    if(p > 3){
      pnames <- c(
        pnames,
        paste0("beta", 1:(p - 3)),
        "phi"
      )
    } else {
      pnames <- c(pnames, "phi")
    }
  }
  if(fam == "poisson" & p > 2){
    pnames <- c(
      pnames,
      paste0("beta", 1:(p - 2))
    )
  }
  
  if(sum(is.na(y)) == 0){
    
    # fit the model
    fit <- optim(init_theta, ricker_count_neg_ll, y = y, fam = fam, X = X)
    theta_hat <- fit$pars
    V <- solve(fit$hessian)
    
    # construct object returns with names
    names(theta_hat) <- pnames
    colnames(V) <- pnames; rownames(V) <- pnames
    
    
    return(list(
      estim = theta_hat,
      V = V
    ))
    
  }
  
  # fit with EM algorithm for data with missing obs
  if(sum(is.na(y)) > 0){
    
    fit <- ricker_EM(y = y, init_theta = init_theta, fam = fam, tol = tol, max_iter = max_iter)
    
    V <- suppEM(y = y, Theta = fit$Theta, V_c = fit$V_c, fam = fam, tol = tol, max_iter = max_iter)
    
    theta_hat <- fit$theta
    names(theta_hat) <- pnames
    colnames(V) <- pnames; rownames(V) <- pnames
    
    return(
      estims = theta_hat,
      V = V
    )
  }
  
}



