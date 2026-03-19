
#' Compute the mode of an empirical distribution
#'
#' @param x Vector of samples from the target distribution.
#'
#' @return The mode of the empirical distribution.
#' 
posterior_mode <- function(x){
  d <- density(x)
  return(d$x[which.max(d$y)])
}





#' Metropolis-Hastings sampling of the posterior with block updates
#'
#' @param dat Data (as a list).
#' @param lp Log-joint probability of the model. This is a function that takes in
#' the current value of \eqn{\boldsymbol \theta} as it's first argument and the
#' data list as its second to compute the log-probability of the model.
#' @param q_rng Proposal random number generator. This function has one argument and takes in
#' the current value of \eqn{\boldsymbol \theta}.
#' @param q_lpdf Log proposal density. This function takes in the proposal for \eqn{\boldsymbol \theta}
#' as its first argument and the current value as its second.
#' @param burnin Number of warm-up samples.
#' @param iter Number of samples (iterations) to keep.
#' @param nthin Thinning ratio for the kept samples.
#'
#' @return A list with a matrix called \code{theta} (one row per sample and one column per parameter).
#' The second element in the list is a vector of whether the proposal was accepted in each iteration.
#' 
MH_block_sample <- function(dat, lp, burnin, iter, nthin = 1){
  
  ## ---- Initializing starting values ----
  fit_init <- fit_ricker_cc(dat$y, fam = dat$fam, off_patch = TRUE)
  theta_init <- fit_init$estim
  # convert alpha to log-alpha
  # but take care that alpha was estimated as positive
  if(theta_init["alpha"] < 0){
    if(fit_init$upper[which(names(fit_init$estim) == "alpha")] < 0){
      theta_init["alpha"] <- log(0.001)
    } else{
      theta_init["alpha"] <- log(fit_init$upper[which(names(fit_init$estim) == "alpha")])
    }
  } else{
    theta_init["alpha"] <- log(fit_init$estim["alpha"])
  }
  if(is.infinite(theta_init["psi"])){
    theta_init["psi"] <- log(100)
  } else {
    theta_init["psi"] <- log(theta_init["psi"])
  }
  names(theta_init)[which(names(theta_init) == "alpha")] <- "lalpha"
  names(theta_init)[which(names(theta_init) == "psi")] <- "lpsi"
  
  ## ---- Defining the proposal distribution ----
  # # define negative log-likelihood in order to get Hessian
  # nll <- function(x, y, fam, psi = NULL){
  #   y_t <- c(y[1, "ytm1"], as.double(y[,"yt"]))
  #   X <- cbind(1, y_t)
  #   x["lalpha"] <- exp(x["lalpha"])
  #   if(fam == "neg_binom"){
  #     theta <- c(x, psi)
  #   } else {
  #     theta <- x
  #   }
  #   ricker_count_neg_ll(theta = theta, y = y_t, X = X, fam = fam)
  # }
  # 
  # optim(theta_init[1:(p - 1)], nll, y = dat$y, fam = fam, psi = exp(theta_init[p]), hessian = T, method = "BFGS")
  # 
  # # derive hessian around optimal values
  # hess <- rootSolve::hessian(nll, theta_init[1:2], y = dat$y, fam = "neg_binom", psi = exp(theta_init["lpsi"]))

  # get covariance matrix for proposal distribution based on 
  # hessian
  Sigma <- diag(fit_init$se)
  # use delta method for se of log(alpha)
  alpha_index <- which(names(fit_init$estim) == "alpha")
  # these are just values I found work well for these types of data
  Sigma[alpha_index, alpha_index] <- fit_init$se[alpha_index] * 10
  if(dat$fam == "neg_binom"){
    Sigma[nrow(Sigma), ncol(Sigma)] <- 0.5
  }
  
  # define proposal distribution
  q_rng <- function(theta, Sigma){
    if(length(theta) == 1){
      prop <- rnorm(1, theta, Sigma)
      names(prop) <- names(theta)
      return(prop)
    } else {
      theta_prop <- mvtnorm::rmvnorm(1, theta, Sigma)[1, ]
      names(theta_prop) <- names(theta)
      return(theta_prop)
    }
  }
  
  # define proposal density
  q_lpdf <- function(prop, theta, Sigma){
    mvtnorm::dmvnorm(prop, theta, Sigma, log = T)
  }
  
  
  ## ---- Main loop ----
  S <- (burnin + iter) * nthin
  p <- length(theta_init)
  theta_samps <- matrix(nrow = S, ncol = p)
  colnames(theta_samps) <- names(theta_init)
  theta_samps[1, ] <- theta_init
  lp_samps <- vector(length = S, mode = "double")
  lp_samps[1] <- lp(theta_init, dat)
  accept <- vector(mode = "double", S)
  
  # if adaptive step sizes are turned on
  adjust_at <- c(1, floor(burnin / adapt) * c(1:adapt))
  adjust_counter <- 2
  
  ### ---- block updates for poisson ----
  if(fam == "poisson"){
    for(s in 2:S){
      
      # define components of the MH algorithm
      theta_s <- theta_samps[s - 1, ]
      theta_prop <- q_rng(theta_s, Sigma)
      lp_prop <- lp(theta_prop, dat)
      lp_curr <- lp_samps[s - 1]
      
      # compute the ratio
      mh_ratio <- exp((lp_prop + q_lpdf(theta_s, theta_prop, Sigma)) - (lp_curr + q_lpdf(theta_prop, theta_s, Sigma)))
      
      if(is.nan(mh_ratio)){
        mh_ratio <- 0
      }
      # accept or reject the update
      A <- min(1, mh_ratio)
      accept[s] <- rbinom(1, 1, prob = A)
      if(accept[s] == 1){
        theta_samps[s, ] <- theta_prop
        lp_samps[s] <- lp_prop
      } else{
        theta_samps[s, ] <- theta_s
        lp_samps[s] <- lp_curr
      }
      
    } 
  }
  
  ### ---- Alternating updates for negative binomial ----
  if(fam == "neg_binom"){
    accept <- cbind(accept, accept)
    for(s in 2:S){
      
      # define components of the MH algorithm
      theta_s <- theta_samps[s - 1, ]
      theta_prop_b <- c(q_rng(theta_s[-p], Sigma[-p, -p]), theta_s[p])
      lp_prop_b <- lp(theta_prop_b, dat)
      lp_curr_b <- lp_samps[s - 1]
      
      # compute the ratio
      mh_ratio_1 <- exp((lp_prop_b + q_lpdf(theta_s, theta_prop_b, Sigma)) - (lp_curr_b + q_lpdf(theta_prop_b, theta_s, Sigma)))
      
      if(is.nan(mh_ratio_1)){
        mh_ratio_1 <- 0
      }
      # accept or reject the update
      A_1 <- min(1, mh_ratio_1)
      
      accept[s, 1] <- rbinom(1, 1, prob = A_1)
      if(accept[s, 1] == 1){
        theta_s <- theta_prop_b
        lp_curr_b <- lp_prop_b
      }
      
      # second round of proposal
      theta_prop_psi <- c(theta_s[1:(p - 1)], q_rng(theta_s[p], Sigma[p, p]))
      lp_prop_psi <- lp(theta_prop_psi, dat)
      
      # compute the ratio
      mh_ratio_2 <- exp((lp_prop_psi + q_lpdf(theta_s, theta_prop_psi, Sigma)) - (lp_curr_b + q_lpdf(theta_prop_psi, theta_s, Sigma)))
      
      if(is.nan(mh_ratio_2)){
        mh_ratio_2 <- 0
      }
      # accept or reject the update
      A_2 <- min(1, mh_ratio_2)
      
      accept[s, 2] <- rbinom(1, 1, prob = A_2)
      if(accept[s, 2] == 1){
        theta_samps[s, ] <- theta_prop_psi
        lp_samps[s] <- lp_prop_psi
      } else{
        theta_samps[s, ] <- theta_s
        lp_samps[s] <- lp_curr_b
      }
      
      # if adaptive step size
      if(adapt > 0 & s == adjust_at[adjust_counter]){
        pa_1 <- mean(accept[adjust_at[adjust_counter - 1]:s, 1])
        pa_2 <- mean(accept[adjust_at[adjust_counter - 1]:s, 2])
        if(pa_1 < 0.2){
          sig_new_1 <- diag(Sigma)[1:2] * 0.75 
          diag(Sigma)[1:(p - 1)] <- sig_new_1
        } else if(pa_1 > 0.6){
          sig_new_1 <- diag(Sigma)[1:2] * 1.5
          diag(Sigma)[1:(p - 1)] <- sig_new_1
        }
        
        # adjustments to stepsize for dispersion
        if(pa_2 < 0.2){
          sig_new21 <- diag(Sigma)[p] * 0.75 
          diag(Sigma)[p] <- sig_new_2
        } else if(pa_2 > 0.6){
          sig_new_2 <- diag(Sigma)[p] * 1.5
          diag(Sigma)[p] <- sig_new_2
        }
        # increment counter
        adjust_counter <- adjust_counter + 1
      }
      
    }
  }
  
  # thin out, then return the post-burnin samples
  samps2keep <- seq(1, S, by = nthin)
  theta_samps2keep <- theta_samps[samps2keep, ]
  colnames(theta_samps2keep) <- names(theta_init)
  return(
    list(
      theta = theta_samps2keep[(burnin + 1):(burnin + iter), ],
      acceptance = accept
    )
  )
  
}



step_out <- function(k, theta_curr, pars_afs, h, lp, dat, max_steps_out = 10000){
  
  # random interval placement
  width <- pars_afs$width[k]
  int <- double(length = 2)
  x_low <- runif(1, min = -width, max = 0)
  x_high <- x_low + width
  Gam_k <- pars_afs$Gamma[, k]
  
  # do the lower bound
  h_low <- lp(x_low * Gam_k + theta_curr, dat)
  i <- 0
  n_expand <- 0
  while (i < max_steps_out & h_low > h) {
    x_low <- x_low - width
    h_low <- lp(x_low * Gam_k + theta_curr, dat)
    i <- i + 1
  }
  int[1] <- x_low
  n_expand <- n_expand + i
  
  # now expand upper bound
  h_high <- lp(x_high * Gam_k + theta_curr, dat)
  i <- 0
  while(i < max_steps_out & h_high > h) {
    x_high <- x_high + width
    h_high <- lp(x_high * Gam_k + theta_curr, dat)
    i <- i + 1
  }
  int[2] <- x_high
  n_expand <- n_expand + i
  
  if(n_expand == 0){ n_expand <- 1 }
  
  return(list(int = int, n_expand = n_expand))
  
}



slice_proposals <- function(k, A_curr, theta_curr, pars_afs, h, lp, dat, adapt_A = T){
  
  h_prop <- -Inf
  prop_counter <- 0
  contraction_counter <- 0
  while(h_prop < h){
    q_prop <- runif(1, min = A_curr[1], max = A_curr[2])
    h_prop <- lp(q_prop * pars_afs$Gamma[, k] + theta_curr, dat)
    if(adapt_A & h_prop < h){
      if(q_prop < 0){
        A_curr[1] <- q_prop
      } else {
        A_curr[2] <- q_prop
      }
      contraction_counter <- contraction_counter + 1
    }
    prop_counter <- prop_counter + 1
  }
  if(prop_counter >= time_out){
    stop(paste0("Failed to find suitable proposal for theta_", k))
  }
  theta_prop <- q_prop * pars_afs$Gamma[, k] + theta_curr
  
  return(
    list(theta = theta_prop, contraction_count = contraction_counter)
  )
  
}



auto_tune_afs <- function(theta_init, dat, lp, burnin = 500, tune_w_after = 1, ratio_target = 0.8, stop_after = 2^10 + 1){
    
    p <- length(theta_init)
    
    # ---- Initialize sampler and conduct tuning ----
    pars_afs <- list(
      Gamma = diag(nrow = p),
      width = runif(p, max = 0.1),
      theta_cov = diag(nrow = p)
    )
    
    # set counters
    iter <- 1
    counter_w <- 1
    # counter_cov <- 1
    convergence <- rep(1, p)
    names(convergence) <- paste("theta", 1:p, sep = "_")
    steps_out <- rep(0, p)
    contractions <- rep(0, p)
    while(iter <= stop_after & any(convergence == 1)){
      
      samp_i <- factor_slice_sampler(theta, dat, lp, pars_afs, iter = 1, control = list(adapt_A = T))
      
      # add to running totals
      steps_out <- samp_i$out_steps[1, ] + steps_out
      contractions <- samp_i$in_steps[1, ] + contractions
      
      # xbar <- xbar + samp_i$theta
      # # note that theta is a row vector here
      # sample_cov <- sample_cov + samp_i$theta %*% t(samp_i$theta)
      
      ## ---- Tuning the widths ----
      for(j in 1:p){
        denom <- steps_out[j] + contractions[j]
        if(denom > 0.0){
          ratio <- steps_out[j] / denom
          
          if(ratio == 0.0){
            ratio <- 1 / denom
          }
          
          multiplier <- (ratio / ratio_target)
          
          # modify width
          pars_afs$width[j] <-
            pars_afs$width[j] * multiplier
          
          if(multiplier > 0.9 & multiplier < 1.1){
            convergence[j] <- 0
          }
        }
      }
      # now reset counter for w and increase number of trials before reset
      counter_w <- 0
      steps_out <- rep(0, p)
      contractions <- rep(0, p)
      iter <- iter + 1
      
    } # end while
    
    ## ---- Estimating factors ----
    theta_burn <- factor_slice_sampler(theta_init, dat, lp, pars_afs, iter = 2000)$theta
    # normalize the covariance estimate
    V <- cov(theta_burn)
    
    # estimate factors
    pars_afs$Gamma <- eigen(V)$vectors
    pars_afs$theta_cov <- V
    
    return(pars_afs)
}



factor_slice_sampler <- function(theta_init, dat, lp, pars_afs, iter = 1000, control = list()){
  
  # set controls for other funs if not set
  alg_control <- list(
    #time_out = 100,
    adapt_A = TRUE,
    max_steps_out = 100
  )
  alg_control <- modifyList(alg_control, control)
  
  # initialize trackers
  p <- length(theta_init)
  theta <- matrix(nrow = iter + 1, ncol = p)
  colnames(theta) <- names(theta_init)
  theta[1, ] <- theta_init
  div_trans <- matrix(0, nrow = iter, ncol = p)
  colnames(div_trans) <- names(theta_init)
  
  if(alg_control$adapt_A){
    out_steps <- matrix(0, nrow = iter, ncol = p)
    in_steps <- matrix(0, nrow = iter, ncol = p)
  }
  
  # alg 2 of Tibbits et al. 2013
  for(i in 1:iter){
    
    # sample random height under the current value
    h_i <- lp(theta[i, ], dat) - rexp(1)

    theta_prop <- theta[i, ]
    for(j in 1:p){
      # approximate interval
      out <- step_out(j, theta_prop, pars_afs, h = h_i, lp = lp, dat = dat)
      A_j <- out$int
      # sample from interval
      samp_i <- slice_proposals(j, A_j, theta_prop, pars_afs, h_i, lp, dat, alg_control$adapt_A)
      
      theta_prop[j] <- samp_i$theta[j]
      
      if(alg_control$adapt_A){
        out_steps[i, j] <- out$n_expand
        in_steps[i, j] <- samp_i$contraction_count
      }
      
    }
    theta[i + 1, ] <- theta_prop
  }
  return(
    list(
      theta = theta[-1, ],
      out_steps = out_steps,
      in_steps = in_steps
    )
  )
}




#' Block Metropolis-Hastings within Gibbs sampling of the posterior for Data Augmentation
#'
#' @param dat Data, as a list.
#' @param fill_rng Gibbs sampler to "fill in" the missing values of \code{y}, conditional on
#' the current value of \eqn{\boldsymbol \theta}.
#' @param lp Function to compute the log-probability of the model, conditional on the current
#' value of \eqn{\boldsymbol \theta} and \eqn{y}.
#' @param q_rng Proposal random number generator. This function has one argument and takes in
#' the current value of \eqn{\boldsymbol \theta}. 
#' @param q_lpdf Log proposal density. This function takes in the proposal for \eqn{\boldsymbol \theta}
#' as its first argument and the current value as its second.
#' @param burnin Number of warm-up samples.
#' @param iter Number of samples (iterations) to keep.
#' @param nthin Thinning ratio for the kept samples.
#'
#' @return A list with a matrix called \code{theta} (one row per sample and one column per parameter). 
#' Column names are inherited from \code{theta_init}. The second element in the list is a matrix \code{y}
#' with one row per iteration and a column for each observation of the response. The values in these columns
#' will be fixed for observed values of y, but will vary for unobserved values. The final element in the 
#' list is a vector of whether the proposal was accepted in each iteration.
#' 
MH_Gibbs_DA <- function(dat, fill_rng, lp, burnin, iter, nthin = 1){
  
  # initialize
  S <- (burnin + iter) * nthin
  
  fit_init <- fit_ricker_cc(dat$y)
  # convert alpha to log-alpha
  # but take care that alpha was estimated as positive
  if(fit_init$estim[2] < 0){
    if(fit_init$upper[2] < 0){
      theta_init <- c(fit_init$estim[1], log(0.001))
    } else{
      theta_init <- c(fit_init$estim[1], log(fit_init$upper[2]))
    }
  } else{
    theta_init <- c(fit_init$estim[1], log(fit_init$estim[2]))
  }
  names(theta_init)[2] <- "lalpha"
  y_full_init <- fill_rng(theta_init, dat)
  nll <- function(x, y){
    
    n <- length(y)
    p <- length(x)
    # compute means
    eta <- vector(mode = "double", length = n)
    eta[1] <- log(y[1])
    for(t in 2:n){
      eta[t] <- log(y[t - 1]) + x[1] - y[t - 1] * exp(x[2])
    }
    
    # return the negative log-likelihood
    return(-sum(
      dpois(x = y[2:n], lambda = exp(eta[2:n]), log = T)
    ))
  }
  
  hess <- optim(theta_init, nll, y = y_full_init, hessian = T)$hessian
  
  # Check eigenvalues
  eigenvalues <- eigen(hess)$values
  if (any(eigenvalues < 1e-8)) {
    warning("Hessian matrix is near-singular. Regularizing...")
    hess <- hess + diag(1e-6, nrow(hess))
  }
  
  Sigma <- solve(hess)
  
  # define proposal distribution
  q_rng <- function(theta, Sigma){
    theta_prop <- mvtnorm::rmvnorm(1, theta, Sigma)[1, ]
    names(theta_prop) <- names(theta)
    return(theta_prop)
  }
  
  # define proposal density
  q_lpdf <- function(prop, theta_s, Sigma){
    mvtnorm::dmvnorm(prop, theta_s, Sigma, log = T)
  }
  
  # preserve names
  theta_samps <- matrix(nrow = S, ncol = length(theta_init))
  colnames(theta_samps) <- names(theta_init)
  theta_samps[1, ] <- theta_init
  
  # tracker for incomplete data
  y_samps <- matrix(nrow = S, ncol = length(dat$y))
  y_samps[1, ] <- y_full_init
  
  # compute initial log probability
  lp_samps <- vector(length = S, mode = "double")
  lp_samps[1] <- lp(theta_init, dat, y_samps[1, ])
  
  # check validity of initial values
  i <- 1
  while(is.infinite(lp_samps[1]) & i < 100){
    theta_init_rt <- runif(length(theta_init), -1, 1)
    names(theta_init_rt) <- names(theta_init)
    y_samps[1, ] <- fill_rng(theta_init_rt, dat)
    lp_samps[1] <- lp(theta_init_rt, dat, y_samps[1, ])
    i <- i + 1
  }
  if(is.infinite(lp_samps[1])){
    stop("Failed to find suitable starting values for the parameters.")
  }
  accept <- vector(mode = "double", S)
  accept2 <- accept
  dat_s <- dat
  
  for(s in 2:S){
    
    # define components of the MH algorithm
    theta_s <- theta_samps[s - 1, ]
    theta_prop <- q_rng(theta_s, Sigma)
    
    dat_s$y <- y_samps[s - 1, ]
    lp_prop <- lp(theta_prop, dat_s)
    lp_curr <- lp_samps[s - 1]
    
    # compute the ratio
    mh_ratio <- exp(
      (lp_prop + q_lpdf(theta_s, theta_prop, Sigma)) - 
        (lp_curr + q_lpdf(theta_prop, theta_s, Sigma))
    )
    
    if(is.nan(mh_ratio)){
      mh_ratio <- 0
    }
    # accept or reject the update
    A <- min(1, mh_ratio)
    accept[s] <- rbinom(1, 1, prob = A)
    if(accept[s] == 1){
      theta_samps[s, ] <- theta_prop
      lp_samps[s] <- lp_prop
    } else{
      theta_samps[s, ] <- theta_s
      lp_samps[s] <- lp_samps[s - 1]
    }
    
    # now propose new y_miss
    y_full_prop <- fill_rng(theta_samps[s, ], dat)
    lp_prop2 <- compute_lp(theta_samps[s, ], dat, y_full = y_full_prop)
    lp_curr2 <- lp_samps[s]
    mh_ratio2 <- exp((lp_prop2) - (lp_curr2))
    
    # accept or reject the update
    A2 <- min(1, mh_ratio2)
    accept2[s] <- rbinom(1, 1, prob = A2)
    if(accept2[s] == 1){
      y_samps[s, ] <- y_full_prop
      lp_samps[s] <- lp_prop2
    } else{
      y_samps[s, ] <- y_samps[s - 1, ]
    }
    
  }
  
  # thin out, then return the post-burnin samples
  samps2keep <- seq(1, S, by = nthin)
  theta_samps2keep <- theta_samps[samps2keep, ]
  y_samps2keep <- y_samps[samps2keep, ]
  return(
    list(
      theta = theta_samps2keep[(burnin + 1):(burnin + iter), ],
      y = y_samps2keep[(burnin + 1):(burnin + iter), ],
      acceptance = accept
    )
  )
  
}





#' Fit a Bayesian Ricker population growth model with (potentially) missing data
#'
#' @param y A vector of counts with NA's in place of any missing data
#' @param fam Character of either \code(c("poisson", "neg_binom"))
#' @param chains Number of MCMC chains to run in parallel. Each will generate \code{samples}
#' samples from the posterior.
#' @param samples Number of posterior samples to keep from each chain
#' @param burnin Number of samples to use as warmup.
#' @param priors_list Named list defining the priors for \eqn{r} and \eqn{\alpha}
#' @param nthin Number of steps between samples that actually get returned.
#' @param return_y Logical to indicate whether to return the posterior samples of the missing
#' observations. This will return posterior modes and credible intervals for all observations,
#' so the observed observations will simply be what was observed with no variance.
#'
#' @return A list of posterior estimates, sds, and CredIs
#'
fit_ricker_DA <- function(
    y, fam = "poisson", 
    chains = 4, 
    samples = 1000, 
    burnin = 2000,
    priors_list = list(
      m_r = 0,
      sd_r = 2.5,
      m_lalpha = -3,
      sd_lalpha = 1,
      m_lpsi = 0,
      sd_lpsi = 100
    ),
    nthin = 5, 
    return_y = FALSE,
    off_patch = FALSE
){
  require(parallel)

  # Check for population extinction
  if(sum(y==0,na.rm=T)>0){
    warning("population extinction caused a divide by zero problem, returning NA")
    return(list(
      NA,
      cause = "population extinction"
    ))
  }
  
  # Check for NaN
  if(any(is.nan(unlist(y)),na.rm=T)){
    warning("NaN found, recode missing data as NA, returning NA")
    return(list(
      NA,
      reason = "NaN found"
    ))
  }
  
  # Check for Inf
  if(any(is.infinite(unlist(y)),na.rm=T)){
    warning("infinite population detected, recheck data returning NA")
    return(list(
      NA,
      reason = "population explosion"
    ))
  }
  
  if(isFALSE(off_patch)){
    # remove starting NAs
    if(is.na(y[1])){
      warning("Removing starting NAs...")
      start <- min(which(!is.na(y)))
      y <- y[start:length(y)]
    } 
  }
  
  # check for too few non-missing sets y(t) and y(t-1) for initial estimates from fit_ricker_cc 
  # compile into sliced dataframe
  if(off_patch){
    n <- nrow(y) + 1
    dat <- y[, c("yt", "ytm1")]
  } else {
    n <- length(y)
    dat <- data.frame(
      yt = y[2:n],
      ytm1 = y[1:(n - 1)]
    )
  }
  
  # check that there are enough complete cases
  ncc <- nrow(dat[complete.cases(dat), ])
  
  if(ncc < 5){
    warning("There are not enough non-missing sets y(t) and y(t-1)")
    return(list(
      NA,
      reason = "missingness limit"
    ))
  }
  
  # ---- Initializing starting values ----
  fit_init <- fit_ricker_cc(dat, fam = fam, off_patch = TRUE)
  theta_init <- fit_init$estim
  # convert alpha to log-alpha
  # but take care that alpha was estimated as positive
  if(theta_init["alpha"] < 0){
    if(fit_init$upper[which(names(fit_init$estim) == "alpha")] < 0){
      theta_init["alpha"] <- log(0.001)
    } else{
      theta_init["alpha"] <- log(fit_init$upper[which(names(fit_init$estim) == "alpha")])
    }
  } else{
    theta_init["alpha"] <- log(fit_init$estim["alpha"])
  }
  if(is.infinite(theta_init["psi"])){
    theta_init["psi"] <- log(100)
  } else {
    theta_init["psi"] <- log(theta_init["psi"])
  }
  names(theta_init)[which(names(theta_init) == "alpha")] <- "lalpha"
  names(theta_init)[which(names(theta_init) == "psi")] <- "lpsi"
  
  
  # create internal function to compute the log-probability
  compute_lp <- function(theta, datlist, y_full = NULL){
    if(is.null(y_full)){
      y_full <- c(datlist$y[1, "ytm1"], as.double(datlist$y[, "yt"]))
    }
    r <- theta["r"]
    alpha <- -exp(theta["lalpha"])
    Xmm <- cbind(
      rep(1, length(y_full)),
      y_full
    )
    if(datlist$fam == "poisson"){
      theta2 = c(r, alpha)
    }
    if(datlist$fam == "neg_binom"){
      psi <- exp(theta["lpsi"])
      theta2 = c(r, alpha, psi = psi)
    }
    lprob <- -ricker_count_neg_ll(theta = theta2, y = y_full, X = Xmm, fam = datlist$fam) +
      dnorm(theta["r"], mean = datlist$m_r, sd = datlist$sd_r, log = T) + 
      dnorm(theta["lalpha"], mean = datlist$m_lalpha, sd = datlist$sd_lalpha, log = T)
    
    if(datlist$fam == "neg_binom"){
      lprob <- lprob +
        dnorm(theta["lpsi"], mean = datlist$m_lpsi, sd = datlist$sd_lpsi, log = T)
    }
    return(unname(lprob))
  }
  
  # function to "fill in" the missing data based on current values of params
  fill_rng <- function(theta, datlist, family = fam){
    y <- datlist$y
    for(t in 1:nrow(y)){
      if(is.na(y[t, "yt"])){
        mu_t <- ricker_step(
          theta = c(theta["r"], -exp(theta["lalpha"])), 
          Nt = y[t, "ytm1"]
        )
        if(family == "poisson"){
          y[t, "yt"] <- zt_poisson_rng(1, mu_t)
        }
        if(family == "neg_binom"){
          y[t, "yt"] <- zt_neg_binom_rng(1, size = exp(theta["lpsi"]), mu = mu_t)
        }
      }
    }
    return(y)
  }
  
  # proposal distribution is now determined based
  
  # compile data
  datlist <- c(
    list(
      y = dat,
      fam = fam
    ),
    priors_list
  )
  
  # Starting point now determined based on empirical Bayes approach
  
  prop_miss <- mean(is.na(dat[, "yt"]))

  if(prop_miss == 0){
    force(ls(envir = environment()))
    cl <- parallel::makeCluster(chains)
    parallel::clusterEvalQ(
      cl,
      {
        source(here::here("Functions/ricker_drop_function.R"))
        source(here::here("Functions/ricker_count_MCMC.R"))
        source(here::here("Functions/ricker_count_likelihood_functions.R"))
        }
    )

    # export the remaining variables

    parallel::clusterExport(
      cl, 
      varlist = ls(envir = environment()), 
      envir = environment()
    )
    
    post_samps <- parallel::clusterCall(
      cl,
      MH_block_sample,
      dat = dat,
      lp = compute_lp,
      burnin = burnin,
      iter = samples,
      nthin = nthin
    )
    parallel::stopCluster(cl)
  }
  
  if(prop_miss > 0){
    force(ls(envir = environment()))
    cl <- parallel::makeCluster(chains)
    parallel::clusterEvalQ(
      cl,
      {
        source(here::here("Functions/ricker_drop_function.R"))
        source(here::here("Functions/ricker_count_MCMC.R"))
        source(here::here("Functions/ricker_count_likelihood_functions.R"))
        }
    )
    
    # export the remaining variables
    parallel::clusterExport(
      cl, 
      varlist = ls(envir = environment()), 
      envir = environment()
    )
    
    post_samps <- parallel::clusterCall(
      cl,
      MH_Gibbs_DA,
      dat = dat,
      lp = compute_lp,
      fill_rng = fill_rng,
      burnin = burnin,
      iter = samples,
      nthin = nthin
    )
    parallel::stopCluster(cl)
    
  } 
  
  # compute rhat statistic for samples
  post_r <- Reduce(
    cbind,
    lapply(post_samps, function(x){x$theta[,"r"]})
  )
  post_lalpha <- Reduce(
    cbind,
    lapply(post_samps, function(x){x$theta[,"lalpha"]})
  )
  
  ses <- apply(cbind(post_r, post_lalpha), 2, sd)
  if(any(ses == 0)){
    stuck_r <- which(ses[1:chains] == 0)
    stuck_alpha <- which(ses[(chains + 1):(2 * chains)] == 0)
    mess <- paste(
      "MCMC sampler got stuck. Chain", 
      stuck_r, "for param r, and chain", 
      stuck_alpha, "for param alpha.\n"
    )
    warning(mess)
    return(list(
      NA,
      reason = "Sampler got stuck"
    ))
  }
  
  # combine posterior samps of r and alpha
  theta_samps <- Reduce(
    rbind,
    lapply(post_samps, function(x){
      x$theta
    })
  )
  theta_samps[,"lalpha"] <- exp(theta_samps[,"lalpha"])
  colnames(theta_samps)[1:2] <- c("r", "alpha")
  if(fam == "neg_binom"){
    theta_samps["lpsi"] <- exp(theta_samps["lpsi"])
    colnames(theta_samps)[3] <- "psi"
  }
  
  # if we want to return posterior estims for missing obs
  if(isTRUE(return_y) & prop_miss > 0){
    y_samps <- Reduce(
      rbind,
      lapply(post_samps, function(x){
        x$y
      })
    )
    theta_samps <- cbind(
      theta_samps, y_samps
    )
    colnames(theta_samps)[3:ncol(theta_samps)] <- 
      paste0("y", 1:n)
  }
  
  # return summaries
  return(
    list(
      estim = apply(theta_samps, 2, mean),
      rhat = c(posterior::rhat(post_r), posterior::rhat(post_lalpha)),
      se = apply(theta_samps, 2, sd),
      lower = apply(theta_samps, 2, quantile, probs = 0.025),
      upper = apply(theta_samps, 2, quantile, probs = 0.975)
    )
  )
  
}

