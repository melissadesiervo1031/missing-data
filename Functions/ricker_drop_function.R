

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
      estim = c(estims, psi = exp(lpsi_curr)),
      se = c(ses, psi = NA_real_),
      lower = c(as.double(cis[, 1]), psi = NA_real_),
      upper = c(as.double(cis[, 2]), psi = NA_real_),
      convergence = ifelse(i < max_iter, 0, 1)
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
    fit <- optim(init, ricker_count_neg_ll_cnstr, y = dat_cc, fam = "poisson", hessian = T)
    
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
  }
  
  # ---- or fit with negbinom ----
  if(fam == "neg_binom"){
    init <- c(r = 1, lalpha = log(0.01), lpsi = log(4))
    fit <- ricker_count_nb_fit(init, dat_cc)
  }
  

  # no need for 1000 betas for projection CI
  if(pro_conf=="none"){
    # return as a list
    return(list(
      estim = estims,
      se = ses,
      lower = as.double(cis[, 1]),
      upper = as.double(cis[, 2]),
      convergence = fit$convergence
    ))
  }
  
  # simulated 1000 betas for projection CI
  if(pro_conf=="sim"){
    sim_per_fit=1000
    
    # draws from Fiducial distribution
    CI_results=mvrnorm(sim_per_fit, fit$par, V)  # Generate random samples
    if("lalpha" %in% colnames(CI_results)){
      CI_results[, "lalpha"] <- exp(CI_results[, "lalpha"])
      colnames(CI_results) <- names(estims)
    }
    
    # return as a list
    return(list(
      estim = estims,
      se = ses,
      lower = as.double(cis[, 1]),
      upper = as.double(cis[, 2]),
      CI_results=CI_results
    ))
    
  }
  
  # bootstrapped 1000 betas for projection CI
  if(pro_conf=="boot"){
    sim_per_fit=1000
    
    boot = function(x, y, f){
      data_sample = y[sample(seq_len(nrow(y)), replace = TRUE), ]
      init <- c(r = 1, lalpha = log(0.01))
      fit_s <- optim(init, ricker_count_neg_ll_cnstr, y = data_sample, fam = f, hessian = T)
      return(fit_s$par)
    }
    
    CI_results <- do.call(rbind, lapply(1:sim_per_fit, boot, dat_cc, f = fam))
    CI_results[, "lalpha"] <- exp(CI_results[, "lalpha"])
    colnames(CI_results) <- names(estims)
    
    # return as a list
    return(list(
      estim = estims,
      se = ses,
      lower = as.double(cis[, 1]),
      upper = as.double(cis[, 2]),
      CI_results=CI_results
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
fit_ricker_drop <- function(y, fam = "poisson", pro_conf="none", off_patch=F){
  
  if(off_patch){ # we have taken in already offset patch data
    y0=y
    y=as.numeric(y0[,1])
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
    
    patches=unique(y0$patchv)
    dat=data.frame(matrix(data=NA,nrow=0,ncol=2))
    for(i in 1:length(patches)){
      patch_i=which(y0$patchv==patches[i])
      print(patch_i)
      y_cc=y0$ytm1[patch_i][complete.cases(y0$ytm1[patch_i])]
      print(y_cc)
      n=length(y_cc)
      if(n<2){
        # skip
      } else {
        dat_i <- data.frame(
          yt = y_cc[2:n],
          ytm1 = y_cc[1:(n - 1)]
        )
        print(dat_i)
        dat=rbind(dat,dat_i)
      }

    }
    
    colnames(dat)=c("yt","ytm1")
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
    fit <- glm(
      yt ~ ytm1, data = dat, 
      family = poisson,
      offset = log(ytm1)
    )
  }
  
  # or fit with negbinom
  if(fam == "neg_binom"){
    fit <- MASS::glm.nb(
      yt ~ ytm1 + offset(log(ytm1)), 
      data = dat
    )
  }
  
  # compile objects and rename
  estims <- coef(fit)
  names(estims) <- c("r", "alpha")
  cis <- confint(fit)
  rownames(cis) <- names(estims)
  ses <- sqrt(diag(vcov(fit)))
  names(ses) <- names(estims)
  

  
  # no need for 1000 betas for projection CI
  if(pro_conf=="none"){
    # return as a list
    return(list(
      estim = estims * c(1, -1),
      se = ses,
      lower = as.double(diag(cis) * c(1, -1)),
      upper = as.double(
        c(cis[1,2], cis[2,1]) * c(1, -1)
      )
    ))
  }
  
  # simulated 1000 betas for projection CI
  if(pro_conf=="sim"){
    sim_per_fit=1000
    
    betas <- coef(fit)     # Get coefficients
    vcov_mat <- vcov(fit)  # Get variance-covariance matrix
    CI_results=mvrnorm(sim_per_fit, betas, vcov_mat)  # Generate random samples
    
    
    # return as a list
    return(list(
      estim = estims * c(1, -1),
      se = ses,
      lower = as.double(diag(cis) * c(1, -1)),
      upper = as.double(
        c(cis[1,2], cis[2,1]) * c(1, -1)
      ),
      CI_results=CI_results
    ))
    
  }
  
  # bootstrapped 1000 betas for projection CI
  if(pro_conf=="boot"){
    sim_per_fit=1000
    
    boot = function(x, model){
      data = model.frame(model)
      data_sample = data[sample(seq_len(nrow(data)), replace = TRUE), ]
      names(data_sample)=c("yt","ytm1","offset")
      coef(update(model, data = data_sample))
    }
    
    CI_results <- do.call(rbind, lapply(1:sim_per_fit, boot, fit))
   
    
    # return as a list
    return(list(
      estim = estims * c(1, -1),
      se = ses,
      lower = as.double(diag(cis) * c(1, -1)),
      upper = as.double(
        c(cis[1,2], cis[2,1]) * c(1, -1)
      ),
      CI_results=CI_results
    ))
  }
  
}





