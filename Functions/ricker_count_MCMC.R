


#' Generate truncated-normal random variates
#'
#' @param n Number of samples to draw
#' @param mu Location parameter for the normal distribution
#' @param sd Standard deviation of the half-normal distribution
#' @param a Lower bound
#' @param b Upper bound
#'
#' @return Vector of n half-normal random variates
#' 
trunc_norm_rng <- function(n = 1, mu = 0, sd = 1, a = -Inf, b = Inf){
  
  l <- pnorm((a - mu) / sd)
  u <- pnorm((b - mu) / sd)
  
  q <- runif(n = n, min = l, max = u)
  return(
    qnorm(q) * sd
  )
  
}






#' Log-pdf of the truncated-normal distribution
#'
#' @param x 
#' @param a Lower bound
#' @param b Upper bound
#' @param sigma Standard deviation of the half-normal
#'
#' @return Real scalar or vector
#' 
trunc_norm_lpdf <- function(x, mu = 0, sd = 1, a = -Inf, b = Inf){
  
  if(x < a | x > b){
    return(log(0))
  } else{
    return(
      -log(sd) + dnorm((x - mu) / sd, log = T) - log(
        pnorm((b - mu) / sd) - pnorm((a - mu) / sd)
      )
    )
  }
  
}






#' Metropolis-Hastings sampling of the posterior with block updates
#'
#' @param theta_init Named vector of initial values of the parameter vector.
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
MH_block_sample <- function(theta_init, dat, lp, q_rng, q_lpdf, burnin, iter, nthin = 1){
  
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
    theta_prop <- q_rng(theta_s)
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
#' @param theta_init Initial values for parameter vector \eqn{\theta}. 
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
#' Column names are inhereted from \code{theta_init}. The second element in the list is a matrix \code{y}
#' with one row per iteration and a column for each observation of the response. The values in these columns
#' will be fixed for observed values of y, but will vary for unobserved values. The final element in the 
#' list is a vector of whether the proposal was accepted in each iteration.
#' 
MH_Gibbs_DA <- function(theta_init, dat, fill_rng, lp, q_rng, q_lpdf, burnin, iter, nthin = 1){
  
  # initialize
  S <- (burnin + iter) * nthin
  theta_samps <- matrix(nrow = S, ncol = length(theta_init))
  
  # preserve names
  colnames(theta_samps) <- names(theta_init)
  theta_samps[1, ] <- theta_init
  lp_samps <- vector(length = S, mode = "double")
  lp_samps[1] <- lp(theta_init, dat)
  accept <- vector(mode = "double", S)
  
  # tracker for incomplete data
  y_samps <- matrix(nrow = S, ncol = length(dat$y))
  y_samps[1, ] <- fill_rng(theta_init, dat)
  
  for(s in 2:S){
    
    # define components of the MH algorithm
    theta_s <- theta_samps[s - 1, ]
    theta_prop <- q_rng(theta_s)
    y_s <- y_samps[s - 1, ]
    lp_prop <- lp(theta_prop, dat, y_s)
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
    # draw a new sample for y
    y_samps[s, ] <- fill_rng(theta_samps[s, ], dat)
    
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



#' Fit a Ricker population growth model with (potentially) missing data
#'
#' @param y A vector of counts with NA's in place of any missing data
#' @param fam Character of either \code(c("poisson", "neg_binom"))
#' @param stan_mod Optional. Can supply a compiled stan model (MUCH FASTER FOR APPLYING
#' TO MANY DATASETS AT ONE) or can be left NULL and the function will compile the
#' model internally based on \code{fam}.
#' @param X Optional covariate matrix (SHOULD NOT INCLUDE Y OR AN INTERCEPT)
#' @param cores Number of cores to use during sampling
#' @param chains Number of HMC chains
#' @param samples Number of posterior samples to keep from each chain
#' @param ... Additional arguments to pass to the \code{control()} for STAN
#'
#' @return A list of posterior estimates, sds, and CredIs
#'
fit_ricker_DA <- function(
    y, fam = "poisson", 
    chains = 4, 
    samples = 1000, burnin = 1000,
    prior_pars, stepsize = NULL
){
  
  if(is.null(stepsize)){
    stepsize <- prior_pars / 5
  }
  priors_list <- split(unname(prior_pars), names(prior_pars))
  n <- length(y)
  S <- samples + burnin
  miss_ids <- which(is.na(y))
  n_miss <- length(miss_ids)
  
  # create internal function to compute the log-probability
  compute_lp <- function(param_list, dat_list, fam = "poisson"){
    
    result <- with(c(param_list, dat_list), {
      r <- theta["r"]
      alpha <- exp(theta["lalpha"])
      lp <- -ricker_count_neg_ll(theta = c(r, alpha), y = y, X = Xmm) +
        trunc_norm_lpdf(r, sd = sd_r, a = 0) + 
        dnorm(theta["lalpha"], mean = m_lalpha, sd = sd_lalpha, log = T)
      lp
    })
    return(result)
  }
  
  # generate parameter objects
  params <- list(
    r = vector(mode = "double", length = S),
    lalpha = vector(mode = "double", length = S)
  )
  
  # create model matrix
  Xmm <- cbind(
    rep(1, n), y
  )
  
  # initial values
  params$r[1] <- trunc_norm_rng(sd = prior_pars["sd_r"], a = 0)
  params$alpha[1] <- rnorm(sd = prior_pars["sd_alpha"], b = 0)
  
  
  ### Gibbs sampler ###
  param_mat <- do.call(cbind, params)
  P <- ncol(param_mat)
  u <- matrix(runif(S * ncol(param_mat)), nrow = S, ncol = ncol(param_mat))
  lp <- vector(mode = "double", length = S)
  
  # first the case with no missing data
  if(n_miss == 0){
    
    dat_list <- list(y = y, Xmm = Xmm)
    
    # initialize log-probability
    lp_current <- compute_lp(
      param_list = c(list(
        theta = param_mat[1, ]
      ), priors_list),
      dat_list = dat_list
    )
    
    
    # step the gibbs sampler forward
    for(s in 2:S){
      
      pvec_s <- param_mat[s - 1, ]
      prop_vec <- pvec_s
      for(p in 1:P){
        
        # proposal steps
        if(p == 1){
          prop_vec[p] <- exp(log(pvec_s[p]) + rnorm(1) * stepsize[p])
        }
        if(p == 2){
          prop_vec[p] <- -exp(log(-pvec_s[p]) + rnorm(1) * stepsize[p])
        }
        
        # log probability of proposal
        lp_sp <- compute_lp(
          param_list = c(list(
            theta = prop_vec 
          ), priors_list),
          dat_list = dat_list
        )
        
        # if log probability of proposal is greater than that of
        # the previous step, accept, otherwise accept with probability
        # equal to the ratio
        if(lp_sp > lp_current){
          lp_current <- lp_sp
          pvec_s[p] <- prop_vec[p]
        } else if(exp(lp_sp - lp_current) > u[s, p]){
          lp_current <- lp_sp
          pvec_s[p] <- prop_vec[p]
        }
      }
      
      param_mat[s, ] <- pvec_s
      lp[s] <- lp_current
      
    }
  } else{ ### Now the case with missing data ###
    
    
    
  }
  
  if(fam == "neg_binom"){
    post_samps <- cbind(
      post_samps,
      as.data.frame(rstan::extract(mfit, pars = c("psi")))
    )
  }
  
  # return summaries
  return(
    list(
      bayes_estim = apply(post_samps, 2, mean),
      post_mode = apply(post_samps, 2, post_mode),
      sd = apply(post_samps, 2, sd),
      l95 = apply(post_samps, 2, quantile, probs = 0.025),
      u95 = apply(post_samps, 2, quantile, probs = 0.975)
    )
  )
  
}