

#' Fit Ricker count model with Poisson errors
#'
#' Estimates Ricker population model parameters under a Poisson observation
#' model via a single call to \code{optim}. The intra-specific competition
#' coefficient \eqn{\alpha} is constrained to be positive through the log
#' parameterization \code{lalpha = log(alpha)}.
#'
#' The Ricker mean function is \eqn{\mu_t = y_{t-1} \exp(r - \alpha y_{t-1})},
#' where counts at time \eqn{t} are assumed to follow a Poisson distribution
#' with mean \eqn{\mu_t}.
#'
#' @param theta_init Named numeric vector of initial parameter values. Must
#'   contain elements named \code{"r"} (intrinsic growth rate) and
#'   \code{"lalpha"} (log of the intra-specific competition coefficient).
#' @param y Either a numeric vector of population counts through time, or a
#'   data frame with columns \code{yt} (count at time \eqn{t}) and
#'   \code{ytm1} (count at time \eqn{t-1}). If a vector is supplied it is
#'   converted internally to the lagged data frame format.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{\code{estim}}{Named numeric vector of point estimates for
#'       \code{r} and \code{alpha} (back-transformed from \code{lalpha}).}
#'     \item{\code{se}}{Named numeric vector of standard errors. The SE for
#'       \code{alpha} is obtained via the delta method.}
#'     \item{\code{lower}}{Numeric vector of lower 95\% confidence limits for
#'       \code{r} and \code{alpha}.}
#'     \item{\code{upper}}{Numeric vector of upper 95\% confidence limits for
#'       \code{r} and \code{alpha}.}
#'     \item{\code{working}}{Named numeric vector of estimates on the working
#'       (log) scale, i.e. \code{r} and \code{lalpha}.}
#'     \item{\code{V}}{Variance-covariance matrix of the working-scale
#'       estimates, derived from the inverse Hessian returned by
#'       \code{optim}.}
#'     \item{\code{convergence}}{Integer convergence flag from \code{optim}:
#'       \code{0} indicates successful convergence.}
#'   }
#'
#' @examples
#' y <- readRDS("data/missingDatasets/pois_sim_randMiss_A.rds")[[1]]$y[[1]]
#' init <- c(r = 1, lalpha = log(0.01))
#' ricker_count_pois_fit(theta_init = init, y = y)
#'
ricker_count_pois_fit <- function(theta_init, y){
  
  fit <- optim(theta_init, ricker_count_neg_ll_cnstr, y = y, fam = "poisson", hessian = T)
  
  # compile results
  estims <- fit$par
  estims["lalpha"] <- exp(estims["lalpha"])
  names(estims)[2] <- "alpha"
  
  # CIs and SEs
  V <- solve(fit$hessian)
  ses <- sqrt(diag(V))
  cis <- mapply(FUN = function(x, s){ x + c(-1, 1) * 2 * s }, x = fit$par, s = ses)
  # convert lalpha
  cis[, 2] <- exp(cis[, 2])
  cis <- t(cis)
  colnames(cis) <- c("lower_95", "upper_95")
  
  # need to use the delta method for se of alpha
  ses[2] <- sqrt(sqrt(ses[2]) * exp(fit$par[2])^2)
  names(ses) <- names(estims)
  
  return(
    list(
      estim = estims,
      se = ses,
      lower = as.double(cis[, 1]),
      upper = as.double(cis[, 2]),
      working = fit$par,
      V = V,
      convergence = fit$convergence,
      nll = fit$value
    )
  )
  
}




#' Fit Ricker count model with negative binomial errors via alternating optimization
#'
#' Estimates Ricker population model parameters under a negative binomial
#' observation model using a two-step alternating optimization scheme.
#' In step 1, the intrinsic growth rate (\code{r}) and log-transformed
#' intra-specific competition coefficient (\code{lalpha}) are optimized
#' conditional on the current overdispersion estimate. In step 2, the
#' log-overdispersion parameter (\code{lpsi}) is optimized conditional on
#' the current predicted means. The two steps alternate until parameter
#' estimates converge or the iteration limit is reached.
#'
#' The Ricker mean function is \eqn{\mu_t = y_{t-1} \exp(r - \alpha y_{t-1})},
#' where \eqn{\alpha} is constrained to be positive via the log
#' parameterization \code{lalpha = log(alpha)}.
#'
#' @param theta Named numeric vector of initial parameter values. Must contain
#'   elements named \code{"r"} (intrinsic growth rate), \code{"lalpha"}
#'   (log of the intra-specific competition coefficient), and \code{"lpsi"}
#'   (log of the negative binomial overdispersion parameter).
#' @param y Either a numeric vector of population counts through time, or a
#'   data frame with columns \code{yt} (count at time \eqn{t}) and
#'   \code{ytm1} (count at time \eqn{t-1}). If a vector is supplied it is
#'   automatically converted to the lagged data frame format.
#' @param tol Numeric scalar. Convergence tolerance; the Euclidean distance
#'   between successive parameter vectors must fall below this value for the
#'   algorithm to be considered converged. Defaults to \code{1e-5}.
#' @param max_iter Integer. Maximum number of alternating-optimization
#'   iterations before the algorithm stops regardless of convergence.
#'   Defaults to \code{500}.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{\code{estim}}{Named numeric vector of point estimates for
#'       \code{r}, \code{alpha} (back-transformed from \code{lalpha}), and
#'       \code{psi} (back-transformed from \code{lpsi}).}
#'     \item{\code{se}}{Named numeric vector of standard errors for \code{r}
#'       and \code{alpha}. The SE for \code{alpha} is obtained via the delta
#'       method. \code{psi} SE is \code{NA} as it is estimated via
#'       \code{optimize()} without a Hessian.}
#'     \item{\code{lower}}{Numeric vector of lower 95\% confidence limits for
#'       \code{r} and \code{alpha}. \code{psi} lower limit is \code{NA}.}
#'     \item{\code{upper}}{Numeric vector of upper 95\% confidence limits for
#'       \code{r} and \code{alpha}. \code{psi} upper limit is \code{NA}.}
#'     \item{\code{convergence}}{Integer flag: \code{0} if the algorithm
#'       converged within \code{max_iter} iterations, \code{1} otherwise.}
#'   }
#'
#' @examples
#' y <- readRDS("data/missingDatasets/nb_sim_randMiss_A.rds")[[1]]$y[[1]]
#' init <- c(r = 1, lalpha = log(0.01), lpsi = log(5))
#' ricker_count_nb_fit(theta = init, y = y)
#'
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
    
    # constrained alpha
    alpha <- exp(theta["lalpha"])
    
    # compute means
    eta <- vector(mode = "double", length = nrow(X) + 1)
    eta[1] <- log(X[1, "ytm1"])
    
    # compute means
    for(t in 1:nrow(X)){
      eta[t + 1] <- log(X[t, "ytm1"]) + theta["r"] - alpha * X[t, "ytm1"]
    }
    
    # return the negative log-likelihood
    return(-sum(
      dnbinom(x = X[["yt"]], mu = exp(eta[2:length(eta)]), size = psi, log = T)
    ))
  }
  
  ## ---- Step 2 ----
  step_2_obj_fun <- function(lpsi, X, mu){
    # return the negative log-likelihood
    return(-sum(
      dnbinom(x = X[["yt"]], mu = mu, size = exp(lpsi), log = T)
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
    eta <- vector(mode = "double", length = nrow(X) + 1)
    eta[1] <- log(X[1, "ytm1"])
    
    # compute means
    for(t in 1:nrow(X)){
      eta[t + 1] <- log(X[t, "ytm1"]) + theta_prop["r"] - exp(theta_prop["lalpha"]) * X[t, "ytm1"]
    }
    
    fit_i2 <- optimize(step_2_obj_fun, interval = c(-5, 5), X = X, mu = exp(eta))
    
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
      estim = c(estims, psi = exp(lpsi_curr)),
      se = c(ses, psi = NA_real_),
      lower = c(as.double(cis[, 1]), psi = NA_real_),
      upper = c(as.double(cis[, 2]), psi = NA_real_),
      working = c(theta_curr, lpsi_curr),
      V = V,
      convergence = ifelse(i < max_iter, 0, 1),
      nll = fit_i1$value
    )
  )
  
}


#' Fit Ricker count model with complete cases only
#'
#' @param y Vector of population counts through time. NA should be used for missing data.
#' @param fam Error distribution to use. Options include c("poisson", "neg_binom").
#' @param pro_conf what type of projection confidence to provide, defaults to "none", "sim" returns 1000 simulated estimates, "boot" returns 1000 bootstrapped estimates
#' 
#' @return List of intrinsic growth factor and intra-specific competitive effect estimates,
#' standard errors, and 95% confidence limits.
#'
#' @examples
#' 
#' y <- readRDS("data/missingDatasets/pois_sim_randMiss_A.rds")[[1]]$y[[1]]
#' fit_ricker_cc(y)
#' 
fit_ricker_cc <- function(y, fam = "poisson", pro_conf="none", off_patch=F){
  
  if(off_patch){ # we have taken in already offset patch data
    y0=y
    y=as.numeric(
      c(y0[1, "ytm1"], y0[,"yt"])
    )
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
  
  
  
  # some useful variables
  n <- length(y)
  
  if(off_patch){ # we have taken in already offset patch data
    dat=y0
  } else {
    # compile into sliced dataframe
    dat <- data.frame(
      yt = y[2:n],
      ytm1 = y[1:(n - 1)]
    )
  }
  
  # drop incomplete cases
  dat_cc <- dat[complete.cases(dat), ]
  
  # check for sufficient complete cases
  if(nrow(dat_cc) < 5){
    warning("There are not enough non-missing sets y(t) and y(t-1)")
    return(list(
      NA,
      reason = "missingness limit"
    ))
  }
  
  # ---- fit with poisson ----
  if(fam == "poisson"){
    init <- c(r = 1, lalpha = log(0.01))
    fit <- ricker_count_pois_fit(init, dat_cc)
  }
  
  # ---- or fit with negbinom ----
  if(fam == "neg_binom"){
    init <- c(r = 1, lalpha = log(0.01), lpsi = log(4))
    fit <- ricker_count_nb_fit(init, dat_cc)
  }
  

  # no need for 1000 betas for projection CI
  if(pro_conf=="none"){
    # return as a list
    return( fit[-which(names(fit) %in% c("working", "V"))] )
  }
  
  # ---- Many betas ----
  ## ---- Fiducial draws ----
  if(pro_conf=="sim"){
    sim_per_fit=1000
    
    # draws from Fiducial distribution
    # put alpha back on log scale
    CI_results=mvrnorm(sim_per_fit, fit$working, fit$V)  # Generate random samples
    if("lalpha" %in% colnames(CI_results)){
      CI_results[, "lalpha"] <- exp(CI_results[, "lalpha"])
      colnames(CI_results) <- names(fit$estim)
    }
    
    # return as a list
    return(c(
      fit[-which(names(fit) %in% c("working", "V"))],
      list(CI_results=CI_results)
    ))
    
  }
  
  ## ---- bootstrapped betas ----
  if(pro_conf=="boot"){
    sim_per_fit=1000
    
    boot = function(x, y, f){
      data_sample = y[sample(seq_len(nrow(y)), replace = TRUE), ]
      init <- c(r = 1, lalpha = log(0.01))
      if(f == "neg_binom"){
        init <- c(init, lpsi = log(4))
        fit_s <- ricker_count_nb_fit(init, data_sample)
      } else {
        fit_s <- ricker_count_pois_fit(init, data_sample)
      }
      return(fit_s$estim)
    }
    
    CI_results <- do.call(rbind, lapply(1:sim_per_fit, boot, dat_cc, f = fam))
    
    # return as a list
    return(c(
      fit[-which(names(fit) %in% c("working", "V"))],
      list(CI_results=CI_results)
    ))
  }
  
}



#' Fit Ricker count model by naively dropping NAs
#' 
#' This function, rather than constructing a complete cases dataframe,
#' uses a naive dropping approach, thus violating the assumption of equal
#' spacing between observations in the time series.
#'
#' @param y Vector of population counts through time. NA should be used for missing data.
#' @param fam Error distribution to use. Options include c("poisson", "neg_binom").
#' @param pro_conf what type of projection confidence to provide, defaults to "none", "sim" returns 1000 simulated estimates, "boot" returns 1000 bootstrapped estimates
#'
#' @return List of intrinsic growth factor and intra-specific competitive effect estimates,
#' standard errors, and 95% confidence limits.
#'
#' @examples
#' 
#' y <- readRDS("data/missingDatasets/pois_sim_randMiss_A.rds")[[1]]$y[[1]]
#' fit_ricker_cc(y)
#' 
fit_ricker_drop <- function(y, fam = "poisson", pro_conf="none", off_patch=F, patch_col = "patch"){
  
  if(off_patch){ # we have taken in already offset patch data
    y0=y
    y=as.numeric(c(y0[1, "ytm1"], y0[, "yt"]))
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
  if(any(is.nan(y), na.rm=T)){
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
  
  
  y <- y[complete.cases(y)]
  
  # fail if trimmed time series is too small 
  if (length(y) <=5) {
    warning("Time series with NAs dropped is too short! Model can't fit well")
    return(list(
      NA,
      reason = "ts too short"
    ))
  }
  
  
  
  if(off_patch){ # we have taken in already offset patch data
    
    patches=unique(y0[, patch_col])
    dat_split <- lapply(patches, function(p, df){ df[df[patch_col] == p, ] }, df = y0)
    dat <- data.frame()
    for(y in dat_split){
      y_i <- c(y[1, "ytm1"], y[, "yt"])
      y_i <- y_i[complete.cases(y_i)]
      df_i <- data.frame(
        yt = y_i[2:length(y_i)],
        ytm1 = y_i[1:(length(y_i) - 1)],
        patch_col = unique(y[[patch_col]])
      )
      dat <- rbind(dat, df_i)
    }

  } else {
    n <- length(y)
    # compile into sliced dataframe
    dat <- data.frame(
      yt = y[2:n],
      ytm1 = y[1:(n - 1)]
    )
  }
  
  # ---- fit with poisson ----
  if(fam == "poisson"){
    init <- c(r = 1, lalpha = log(0.01))
    fit <- ricker_count_pois_fit(init, dat)
  }
  
  # ---- or fit with negbinom ----
  if(fam == "neg_binom"){
    init <- c(r = 1, lalpha = log(0.01), lpsi = log(4))
    fit <- ricker_count_nb_fit(init, dat)
  }
  
  
  # no need for 1000 betas for projection CI
  if(pro_conf=="none"){
    # return as a list
    return( fit[-which(names(fit) %in% c("working", "V"))] )
  }
  
  # ---- Many betas ----
  ## ---- Fiducial draws ----
  if(pro_conf=="sim"){
    sim_per_fit=1000
    
    # draws from Fiducial distribution
    # put alpha back on log scale
    CI_results=mvrnorm(sim_per_fit, fit$working, fit$V)  # Generate random samples
    if("lalpha" %in% colnames(CI_results)){
      CI_results[, "lalpha"] <- exp(CI_results[, "lalpha"])
      colnames(CI_results) <- names(fit$estim)
    }
    
    # return as a list
    return(c(
      fit[-which(names(fit) %in% c("working", "V"))],
      list(CI_results=CI_results)
    ))
    
  }
  
  ## ---- bootstrapped betas ----
  if(pro_conf=="boot"){
    sim_per_fit=1000
    
    boot = function(x, y, f){
      data_sample = y[sample(seq_len(nrow(y)), replace = TRUE), ]
      init <- c(r = 1, lalpha = log(0.01))
      if(f == "neg_binom"){
        init <- c(init, lpsi = log(4))
        fit_s <- ricker_count_nb_fit(init, data_sample)
      } else {
        fit_s <- ricker_count_pois_fit(init, data_sample)
      }
      return(fit_s$estim)
    }
    
    CI_results <- do.call(rbind, lapply(1:sim_per_fit, boot, dat_cc, f = fam))
    
    # return as a list
    return(c(
      fit[-which(names(fit) %in% c("working", "V"))],
      list(CI_results=CI_results)
    ))
  }
  
}





