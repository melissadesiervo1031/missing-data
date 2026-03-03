####################################################
# This is a user-defined function to fit Ricker
# population models to count data, potentially with 
# missing observations, using EM
####################################################

source(here::here("Functions/ricker_count_likelihood_functions.R"))
#' Fitting the Ricker population model to count data
#' 
#' This function uses an Expectation Maximization (EM) approach to fit a Ricker population model
#' with Poisson or Negative-Binomial demographic stochasticity.
#'
#' @param y Vector of population counts with NAs placed in positions in which there are missing
#' observations (if any).
#' @param init_theta Initial values for the parameter vector.
#' @param fam Family of the error distribution. Can be either \code{"poisson"} or \code{"neg_binom"}
#' @param tol Tolerance for convergence.
#' @param max_iter Maximum number of iterations to run the EM algorithm before stopping.
#'
#' @return List with estimates and conditional standard errors for \code{theta}, estimates of the
#' full data vector with latent values filled in (\code{z_hat}), and a convergence code. 0 means 
#' the algorithm converged before being stopped.
#' 
ricker_EM <- function(y, init_theta, fam = "poisson", tol = 1e-5, max_iter = 50, off_patch = FALSE){
  
  if(off_patch){
    n <- nrow(y) + 1
  } else{
    n <- length(y)
    # remove starting NAs
    if(is.na(y[1])){
      obs <- which(!is.na(y))
      y <- y[min(obs):n]
      n <- length(y)
    }
  }
  p <- length(init_theta)
  
  # define initial parameter vector
  Theta <- matrix(init_theta, ncol = length(init_theta), nrow= 1)
  
  # define matrix of latent vectors
  Z <- matrix(data = NA, ncol = n, nrow = 1)
  
  nll <- vector(mode = "double")
  
  # iterate through
  dif <- 1
  s <- 1
  
  while(dif > tol & s <= max_iter){
    
    #initialize latent data
    z_s <- y
    
    # determine parameter sets
    theta <- as.double(Theta[s, ])
    if(fam == "poisson"){
      beta <- theta
    }
    if(fam == "neg_binom"){
      beta <- theta[-p]
      phi <- theta[p]
    }
    
    # fill in missing values with their expected values
    # first the single vector option
    if(isFALSE(off_patch)){
      for(t in 2:n){
        if(is.na(z_s[t])){
          z_s[t] <- round(ricker_step(beta, z_s[t - 1]))
          # if rounds to zero, round up instead
          if(z_s[t] == 0){z_s[t] <- 1}
        }
      }
      X <- cbind(
        rep(1, n),
        z_s
      )
    }
    
    if(off_patch) {
      for(t in 1:nrow(z_s)){
        if(is.na(z_s[t, "yt"])){
          z_s[t, "yt"] <- round(ricker_step(beta, z_s[t, "ytm1"]))
          if(t < nrow(z_s)){
            z_s[t + 1, "ytm1"] <- z_s[t, "yt"]
          }
        }
      }
      # round up if rounding to zero
      z_s[z_s == 0] <- 1
      
      # undo the offsetting for use with ricker_count_neg_ll
      z_s_vec <- c(z_s[1, "ytm1"], z_s[, "yt"])
      X <- cbind(
        rep(1, n),
        z_s_vec
      )
      # overwrite for later steps
      z_s <- z_s_vec
    }

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
    nll <- c(nll, fit_s$value)
    
    # calculate the mean relative difference
    if(s == 1){
      dif <- dif
    } else{
      dif <- abs(
        nll[s] - nll[s - 1]
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
#' @param y Vector of population counts, with NA in the place of missing observations. If \code{off_patch == TRUE},
#' this should be a dataframe with \code{yt} and \code{ytm1} in separate columns, pre-filtered such
#' that \code{ytm1} for a given time series or "patch" does not start with \code{NA}.
#' @param fam Error family. Can be either c("poisson", "neg_binom").
#' @param off_patch Logical. Set to \code{TRUE} when the data come from multiple replicate time series.
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
fit_ricker_EM <- function(y, fam = "poisson", off_patch = FALSE, ...){
  
  args <- list(...)
  
  if(off_patch){
    if(!all(c("yt", "ytm1") %in% colnames(y))){
      stop("If feeding in offset data, ensure that yt and ytm1 are named columns.")
    }
    y <- as.matrix(y[,c("ytm1", "yt")])
  }
  # Check for population extinction
  if(sum(y==0,na.rm=T)>1){
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
  
  if(!exists("init_theta", args)){
    init_theta <- c(0.5, -0.01)
    if(fam == "neg_binom"){
      init_theta <- c(init_theta, 5)
    }
    args$init_theta <- init_theta
  }
  if(!exists("tol", args)){args$tol = 1e-5}
  if(!exists("max_iter", args)){args$max_iter = 50}
  
  args <- c(
    list(
      y = y,
      fam = fam,
      off_patch = off_patch
    ),
    args
  )
  
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
  
  if(fit$convergence == 1){
    warning("EM algorithm did not converge. Consider different starting values or increasing max_iter.")
  }
  
  return(list(
    estim = estims,
    se = NA,
    lower = NA,
    upper = NA
  ))
  
}


# NOTE: DG - I don't think we need these. The above should work for pre-formatted data or a single
# time series. If fitting separately, couldn't we just use apply(fit_ricker_EM, list) to a list of
# datasets?

#' EM algorithm for multiple independent Ricker time series
#' 
#' @param y_list List of population count vectors, NAs for missing observations
#' @param init_theta Initial parameter vector
#' @param fam Error family: "poisson" or "neg_binom"
#' @param tol Convergence tolerance
#' @param max_iter Maximum EM iterations
#' 
# ricker_EM_multi <- function(y_list, init_theta, fam = "poisson", tol = 1e-5, max_iter = 50){
#   
#   K <- length(y_list)
#   p <- length(init_theta)
#   
#   # Remove leading NAs and store cleaned series
#   y_list <- lapply(y_list, function(y){
#     if(all(is.na(y))) return(NULL)
#     if(is.na(y[1])){
#       obs <- which(!is.na(y))
#       y <- y[min(obs):length(y)]
#     }
#     return(y)
#   })
#   
#   # Remove completely missing series
#   n_before <- length(y_list)
#   y_list <- Filter(Negate(is.null), y_list)
#   n_after <- length(y_list)
#   
#   if(n_after == 0) stop("All series are completely missing, cannot fit model.")
#   if(n_after < n_before) warning(sprintf("%d completely missing series removed before fitting.", n_before - n_after))
#   
#   # Initialize parameter history
#   Theta <- matrix(init_theta, nrow = 1)
#   
#   dif <- 1
#   s <- 1
#   
#   # Store the filled-in series at each iteration (list of lists)
#   Z_list_history <- list()
#   
#   while(dif > tol & s <= max_iter){
#     
#     theta <- as.double(Theta[s, ])
#     if(fam == "neg_binom"){
#       beta <- theta[-p]
#       phi  <- theta[p]
#     } else {
#       beta <- theta
#     }
#     
#     # ---- E-step: fill in missing values for each series independently ----
#     z_list_s <- vector("list", K)
#     
#     for(i in seq_len(K)){
#       z_s <- y_list[[i]]
#       n_i <- length(z_s)
#       
#       for(t in 2:n_i){
#         if(is.na(z_s[t])){
#           z_s[t] <- max(1, round(ricker_step(beta, z_s[t - 1])))
#         }
#       }
#       z_list_s[[i]] <- z_s
#     }
#     
#     # Build model matrices for each filled series
#     X_list_s <- lapply(z_list_s, function(z){
#       cbind(1, z)
#     })
#     
#     # ---- M-step: optimize pooled likelihood ----
#     fit_s <- tryCatch({
#       optim(
#         par = theta,
#         fn  = ricker_count_neg_ll_multi,
#         y_list = z_list_s,
#         X_list = X_list_s,
#         fam    = fam,
#         hessian = TRUE
#       )
#     }, error = function(e){
#       message("Error in M-step optim: ", conditionMessage(e))
#       return(NA)
#     })
#     
#     if(is.na(fit_s[1])) return(NA)
#     
#     Theta <- rbind(Theta, fit_s$par)
#     Z_list_history[[s]] <- z_list_s
#     
#     # ---- Check convergence: change in pooled NLL ----
#     if(s == 1){
#       dif <- dif
#     } else {
#       nll_prev <- ricker_count_neg_ll_multi(Theta[s, ],     Z_list_history[[s - 1]], X_list_s, fam)
#       nll_curr <- ricker_count_neg_ll_multi(Theta[s + 1, ], z_list_s,                X_list_s, fam)
#       dif <- abs(nll_prev - nll_curr)
#     }
#     
#     s <- s + 1
#   }
#   
#   theta_star  <- as.double(Theta[nrow(Theta), ])
#   convergence <- as.numeric(s > max_iter)
#   
#   return(list(
#     theta      = theta_star,
#     Theta      = Theta,
#     z_list     = z_list_s,       # final filled-in series
#     convergence = convergence
#   ))
# }



#' Wrapper to fit pooled Ricker EM across multiple time series
#' 
#' @param y_list List of count vectors with NAs for missing observations
#' @param fam Error family: "poisson" or "neg_binom"
#' @param ... Additional args passed to ricker_EM_multi (e.g. init_theta, max_iter)
#' 
# fit_ricker_EM_multi <- function(y_list, fam = "poisson", ...){
#   
#   y_list <- lapply(y_list, function(y){
#     if(all(is.na(y))) return(NULL)
#     if(is.na(y[1])){
#       obs <- which(!is.na(y))
#       y <- y[min(obs):length(y)]
#     }
#     return(y)
#   })
#   
#   # Remove completely missing series
#   n_before <- length(y_list)
#   y_list <- Filter(Negate(is.null), y_list)
#   n_after <- length(y_list)
#   
#   if(n_after == 0) stop("All series are completely missing, cannot fit model.")
#   if(n_after < n_before) warning(sprintf("%d completely missing series removed before fitting.", n_before - n_after))
#   
#   # Input checks across all series
#   for(i in seq_along(y_list)){
#     y <- y_list[[i]]
#     if(sum(y == 0, na.rm = TRUE) > 1){
#       warning(sprintf("Series %d: population extinction detected, returning NA", i))
#       return(list(NA, cause = "population extinction"))
#     }
#     if(any(is.nan(y), na.rm = TRUE)){
#       warning(sprintf("Series %d: NaN found, recode missing data as NA", i))
#       return(list(NA, reason = "NaN found"))
#     }
#     if(any(is.infinite(y), na.rm = TRUE)){
#       warning(sprintf("Series %d: infinite population detected", i))
#       return(list(NA, reason = "population explosion"))
#     }
#   }
#   
#   init_theta <- c(0.5, -0.01)
#   if(fam == "neg_binom") init_theta <- c(init_theta, 10)
#   
#   args <- list(
#     y_list     = y_list,
#     fam        = fam,
#     init_theta = init_theta,
#     tol        = 1e-5,
#     max_iter   = 50
#   )
#   args[names(list(...))] <- list(...)
#   
#   fit <- do.call(ricker_EM_multi, args)
#   
#   if(is.na(fit[1])){
#     return(list(NA, reason = "optimization error"))
#   }
#   
#   if(fam == "neg_binom"){
#     estims <- fit$theta * c(1, -1, 1)
#     names(estims) <- c("r", "alpha", "psi")
#   } else {
#     estims <- fit$theta * c(1, -1)
#     names(estims) <- c("r", "alpha")
#   }
#   
#   return(list(
#     estim      = estims,
#     se         = NA,
#     lower      = NA,
#     upper      = NA,
#     convergence = fit$convergence,
#     z_list     = fit$z_list    # filled-in series, useful for diagnostics
#   ))
# }
