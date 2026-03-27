####################################################-
# This is a user-defined function to fit Ricker
# population models to count data, potentially with 
# missing observations, using EM
####################################################-

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
ricker_EM <- function(y, init_theta, fam = "poisson", tol = 1e-5, max_iter = 50){
  
  if(!all(c("yt", "ytm1") %in% colnames(y))){
    stop("Expecting y as a dataframe with yt and ytm1 as named columns.")
  }
  p <- length(init_theta)
  
  # define initial parameter vector
  Theta <- matrix(init_theta, ncol = length(init_theta), nrow= 1)
  colnames(Theta) <- names(init_theta)
  
  # define matrix of latent vectors
  Z <- matrix(data = NA, ncol = nrow(y), nrow = 1)
  
  nll <- vector(mode = "double")
  
  # iterate through
  dif <- 1
  s <- 1
  
  while(dif > tol & s <= max_iter){
    
    #initialize latent data
    z_s <- y
    
    theta_curr <- Theta[s, ]
    if(fam == "neg_binom"){
      beta_curr <- theta_curr[1:2]
    } else {
      beta_curr <- theta_curr
    }
    beta_curr[2] <- exp(beta_curr[2])
    
    
    # fill in missing values with their expected values
    for(t in 1:nrow(z_s)){
      if(is.na(z_s[t, "yt"])){
        z_s[t, "yt"] <- round(ricker_step(beta_curr, z_s[t, "ytm1"]))
        if(z_s[t, "yt"] == 0){
          z_s[t, "yt"] <- 1
        }
        if(t < nrow(z_s)){
          if(is.na(z_s[t + 1, "ytm1"])){
            z_s[t + 1, "ytm1"] <- z_s[t, "yt"]
          }
        }
      }
    }

    # adding try catch here to avoid error "function cannot be evaluated at initial parameters"
    fit_s <- tryCatch({
      if(fam == "poisson"){
        ricker_count_pois_fit(theta_curr, z_s)
      } else if(fam == "neg_binom"){
        ricker_count_nb_fit(theta_curr, z_s)
      } else {
        stop("fam must be poisson or neg_binom")
      }
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
      fit_s$working
    )
    Z <- rbind(
      Z, z_s[["yt"]]
    )
    nll <- c(nll, fit_s$nll)
    
    # calculate the difference
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
  theta_star <- Theta[nrow(Theta), ]
  theta_star["lalpha"] <- exp(theta_star["lalpha"])
  names(theta_star)[which(names(theta_star) == "lalpha")] <- "alpha"
  if(fam == "neg_binom"){
    theta_star["lpsi"] <- exp(theta_star["lpsi"])
    names(theta_star)[which(names(theta_star) == "lpsi")] <- "psi"
  }
  convergence <- as.numeric(s == max_iter)
  
  return(list(
    theta = theta_star,
    Theta = Theta,
    z_s = z_s,
    convergence = convergence,
    nll = nll[length(nll)]
  ))
  
}



#' Fit a Ricker population model to count data using the EM algorithm
#'
#' Fits the stochastic Ricker model to a single time series of population counts,
#' optionally with missing observations. Missing values are handled via Expectation
#' Maximization (EM): at each E-step, missing counts are replaced by their expected
#' value under the current parameter estimates; at each M-step, parameters are
#' re-estimated by maximum likelihood on the filled-in series. Supports Poisson
#' and Negative Binomial error distributions.
#'
#' @param y Either (1) a numeric vector of population counts with \code{NA} marking
#'   missing observations, or (2) if \code{off_patch = TRUE}, a data frame with columns
#'   \code{yt} (count at time \eqn{t}) and \code{ytm1} (count at time \eqn{t-1}),
#'   already formatted into lag pairs and filtered so that no series begins with
#'   \code{ytm1 = NA}. Leading \code{NA}s in a vector input are stripped with a warning.
#' @param fam Error distribution family. One of \code{"poisson"} (default) or
#'   \code{"neg_binom"}.
#' @param off_patch Logical. Set to \code{TRUE} when \code{y} has already been
#'   formatted into lag pairs (i.e., the \code{yt}/\code{ytm1} data frame format).
#'   Use this when fitting to data from multiple replicate time series that have
#'   been row-bound into a single data frame prior to calling this function.
#'   Default is \code{FALSE}.
#' @param ... Additional arguments passed to the internal \code{ricker_EM} algorithm:
#'   \describe{
#'     \item{\code{init_theta}}{Named numeric vector of starting values. For Poisson,
#'       \code{c(r = 0.5, lalpha = log(0.01))}; for Negative Binomial, also includes
#'       \code{lpsi = log(5)}. Note that \code{alpha} is parameterized on the log scale
#'       internally as \code{lalpha = log(alpha)}, and similarly for the NB dispersion
#'       parameter \code{psi}.}
#'     \item{\code{tol}}{Convergence tolerance on the change in negative log-likelihood
#'       between successive EM iterations. Default \code{1e-5}.}
#'     \item{\code{max_iter}}{Maximum number of EM iterations before stopping.
#'       Default \code{50}.}
#'   }
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{\code{estim}}{Named numeric vector of parameter estimates on the natural
#'       scale: \code{r} (intrinsic growth rate) and \code{alpha} (intraspecific
#'       competition coefficient, positive). For \code{fam = "neg_binom"}, also
#'       includes \code{psi} (NB dispersion parameter).}
#'     \item{\code{se}}{Always \code{NA}. The EM algorithm does not directly yield
#'       standard errors; use \code{fit_ricker_DA} for posterior uncertainty
#'       quantification.}
#'     \item{\code{lower}}{Always \code{NA} (see \code{se}).}
#'     \item{\code{upper}}{Always \code{NA} (see \code{se}).}
#'     \item{\code{convergence}}{Convergence code. \code{0} indicates the algorithm
#'       converged (NLL change fell below \code{tol}) before reaching \code{max_iter};
#'       \code{1} indicates the algorithm was stopped at \code{max_iter} without
#'       converging. In the latter case a warning is issued.}
#'   }
#'   Returns a list with \code{NA} as its first element (and a named \code{cause} or
#'   \code{reason} element) if the series contains zeros (population extinction),
#'   \code{NaN}, \code{Inf}, or if the M-step optimizer fails.
#'
#' @seealso \code{\link{ricker_EM}} for the underlying EM implementation,
#'   \code{\link{fit_ricker_DA}} for Bayesian estimation with posterior uncertainty,
#'   \code{\link{fit_ricker_cc}} for complete-case (listwise deletion) MLE.
#'
#' @examples
#' y <- readRDS("data/missingDatasets/pois_sim_randMiss_A.rds")[[1]]$y[[1]]
#' fit_ricker_EM(y)
#'
#' # Negative Binomial with custom starting values and more iterations
#' fit_ricker_EM(y, fam = "neg_binom", init_theta = c(r = 1, lalpha = log(0.05), lpsi = log(10)),
#'               max_iter = 100)
#'
fit_ricker_EM <- function(y, fam = "poisson", off_patch = FALSE, ...){
  
  args <- list(...)
  
  if(off_patch){
    if(!all(c("yt", "ytm1") %in% colnames(y))){
      stop("If feeding in offset data, ensure that yt and ytm1 are named columns.")
    }
    y <- as.matrix(y[,c("yt", "ytm1")])
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
  if(is.vector(y)){
    if(is.na(y[1])){
      warning("Removing starting NAs...")
      start <- min(which(!is.na(y)))
      y <- y[start:length(y)]
    }
    dat <- data.frame(
      yt = y[2:length(y)],
      ytm1 = y[1:(length(y) - 1)]
    )
  } else {
    dat <- as.data.frame(y[, c("yt", "ytm1")])
  }
  
  if(!exists("init_theta", args)){
    init_theta <- c(r = 0.5, lalpha = log(0.01))
    if(fam == "neg_binom"){
      init_theta <- c(init_theta, lpsi = log(5))
    }
    args$init_theta <- init_theta
  }
  if(!exists("tol", args)){args$tol = 1e-5}
  if(!exists("max_iter", args)){args$max_iter = 50}
  
  args <- c(
    list(
      y = dat,
      fam = fam
    ),
    args
  )
  
  fit <- do.call(ricker_EM, args = args)
  
  # adding condition here to avoid error "function cannot be evaluated at initial parameters"
  if(is.na(fit[1])){
    return(list(NA,reason="we have had an error in evaluating function at initial parameters"))
  }
  
  if(fit$convergence == 1){
    warning("EM algorithm did not converge. Consider different starting values or increasing max_iter.")
  }
  
  return(list(
    estim = fit$theta,
    se = NA,
    lower = NA,
    upper = NA,
    convergence = fit$convergence
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
