
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
  
  fit_init <- fit_ricker_cc(dat$y)
  # convert alpha to log-alpha
  # but take care that alpha was estimated as positive
  if(fit_init$estim[2] < 0){
    if(fit_init$upper[2] < 0){
      theta_init <- c(fit_init$estim[1], 0.001)
    } else{
      theta_init <- c(fit_init$estim[1], log(fit_init$upper[2]))
    }
  }
  names(theta_init)[2] <- "lalpha"
  y_full_init <- fill_rng(theta_init, dat)
  # define negative log-likelihood
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
  
  # derive hessian around optimal values
  hess <- optim(theta_init, nll, y = y_full_init, hessian = T)$hessian
  
  # get covariance matrix for proposal distribution based on 
  # hessian
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
  
  
  S <- (burnin + iter) * nthin
  theta_samps <- matrix(nrow = S, ncol = length(theta_init))
  colnames(theta_samps) <- names(theta_init)
  theta_samps[1, ] <- theta_init
  lp_samps <- vector(length = S, mode = "double")
  lp_samps[1] <- lp(theta_init, dat)
  accept <- vector(mode = "double", S)
  
  for(s in 2:S){
    
    # define components of the MH algorithm
    theta_s <- theta_samps[s - 1, ]
    theta_prop <- q_rng(theta_s, sd = 0.02)
    lp_prop <- lp(theta_prop, dat)
    lp_curr <- lp_samps[s - 1]
    
    # compute the ratio
    mh_ratio <- exp((lp_prop + q_lpdf(theta_s, theta_prop)) - (lp_curr + q_lpdf(theta_prop, theta_s)))
    
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
  
  # thin out, then return the post-burnin samples
  samps2keep <- seq(1, S, by = nthin)
  theta_samps2keep <- theta_samps[samps2keep, ]
  return(
    list(
      theta = theta_samps2keep[(burnin + 1):(burnin + iter), ],
      acceptance = accept
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
      theta_init <- c(fit_init$estim[1], 0.001)
    } else{
      theta_init <- c(fit_init$estim[1], log(fit_init$upper[2]))
    }
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
      sd_lalpha = 1
    ),
    nthin = 5, 
    return_y = FALSE
){
  require(parallel)
  if(fam == "neg_binom"){
    stop(
      "Implementation of Negative Binomial model is still in the works."
    )
  }
  
  # Check for population extinction
  if(sum(y==0,na.rm=T)>0){
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
  
  
  # create internal function to compute the log-probability
  compute_lp <- function(theta, dat, y_full = NULL, family = fam){
    if(is.null(y_full)){
      y_full <- dat$y
    }
    r <- theta["r"]
    alpha <- -exp(theta["lalpha"])
    Xmm <- cbind(
      rep(1, length(y_full)),
      y_full
    )
    lp <- -ricker_count_neg_ll(theta = c(r, alpha), y = y_full, X = Xmm, fam = family) +
      dnorm(r, mean = dat$m_r, sd = dat$sd_r, log = T) + 
      dnorm(theta["lalpha"], mean = dat$m_lalpha, sd = dat$sd_lalpha, log = T)
    return(unname(lp))
  }
  
  # function to "fill in" the missing data based on current values of params
  fill_rng <- function(theta, dat, family = fam){
    y <- dat$y
    n <- length(y)
    y_full <- vector("double", length = n)
    y_full[1] <- y[1]
    for(t in 2:n){
      if(is.na(y[t])){
        mu_t <- ricker_step(
          theta = c(theta["r"], -exp(theta["lalpha"])), 
          Nt = y_full[t-1]
        )
        if(family == "poisson"){
          y_full[t] <- rpois(1, mu_t)
        }
        if(family == "neg_binom"){
          y_full[t] <- rnbinom(1, size = exp(theta["lpsi"]), mu = mu_t)
        }
      } else{
        y_full[t] <- y[t]
      }
    }
    # replace zeros (sort of hacky!!!)
    y_full[y_full == 0] <- 1
    return(y_full)
  }
  
  # proposal distribution is now determined based
  # on empirical bayes approach
  # # define proposal distribution
  # q_rng <- function(theta, sd){
  #   p <- length(theta)
  #   theta_prop <- rnorm(p, mean = theta, sd = sd)
  #   names(theta_prop) <- names(theta)
  #   return(theta_prop)
  # }
  # 
  # # define proposal density
  # q_lpdf <- function(prop, theta_s, sd = stepsize){
  #   sum(dnorm(prop, mean = theta_s, sd = sd, log = T))
  # }
  
  # generate some useful variables
  n <- length(y)
  p <- 2
  
  # compile data
  dat <- c(
    list(
      y = y
    ),
    priors_list
  )
  
  # Starting point now determined based on empirical Bayes approach
  # # initialize
  # theta_init <- c(
  #   runif(1, max = 1),
  #   runif(1, min = -4, max = -1)
  # )
  # names(theta_init) <- c("r", "lalpha")
  
  prop_miss <- mean(is.na(y))

  if(prop_miss == 0 & fam == "poisson"){
    force(ls(envir = environment()))
    cl <- parallel::makeCluster(chains)
    parallel::clusterEvalQ(
      cl,
      {
        source(here::here("Functions/ricker_drop_function.R"))
        source(here::here("Functions/ricker_count_MCMC.R"))
        source(here::here("Functions/ricker_count_cc.R"))
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
  
  if(prop_miss > 0 & fam == "poisson" ){
    force(ls(envir = environment()))
    cl <- parallel::makeCluster(chains)
    parallel::clusterEvalQ(
      cl,
      {
        source(here::here("Functions/ricker_drop_function.R"))
        source(here::here("Functions/ricker_count_MCMC.R"))
        source(here::here("Functions/ricker_count_cc.R"))
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
  theta_samps[,2] <- exp(theta_samps[,2])
  colnames(theta_samps) <- c("r", "alpha")
  
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

