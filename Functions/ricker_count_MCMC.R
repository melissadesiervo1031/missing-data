
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



#' Expand a slice interval along a factor direction using the stepping-out procedure
#'
#' Starting from a randomly placed initial interval of width \code{pars_afs$width[k]},
#' this function expands the lower and upper bounds until both lie below the
#' log-probability threshold \code{h}, implementing the stepping-out procedure of
#' Neal (2003) projected onto factor direction \code{k}.
#'
#' @param k Integer index of the factor direction to step out along.
#' @param theta_curr Numeric vector giving the current parameter values.
#' @param pars_afs List of factor slice sampler parameters. Must contain:
#'   \describe{
#'     \item{\code{width}}{Numeric vector of interval widths, one per factor direction.}
#'     \item{\code{Gamma}}{Matrix whose columns are the factor directions (eigenvectors).}
#'   }
#' @param h Numeric scalar. Log-probability threshold defining the current slice height.
#' @param lp Function to evaluate the log-joint probability. Must accept a parameter
#'   vector as its first argument and the data list as its second.
#' @param dat Data list passed to \code{lp}.
#' @param max_steps_out Maximum number of expansion steps allowed in each direction
#'   before giving up. Defaults to \code{10000}.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{int}}{Numeric vector of length 2 giving the lower and upper bounds
#'       of the expanded interval on the standardised factor scale.}
#'     \item{\code{n_expand}}{Total number of expansion steps taken across both
#'       directions (lower and upper combined); set to 1 if no expansion was needed.}
#'   }
#'
step_out <- function(k, theta_curr, pars_afs, h, lp, dat, max_steps_out = 10000){
  
  # random interval placement
  width <- pars_afs$width[k]
  int <- double(length = 2)
  x_low <- runif(1, min = -width, max = 0)
  x_high <- x_low + width
  Gam_k <- pars_afs$Gamma[, k]
  
  # safe_lp <- function(x) {
  #   val <- lp(x, dat)
  #   if (!is.finite(val)) {
  #     warning(sprintf("step_out: lp returned %s at factor %d — treating as -Inf", val, k))
  #     return(-Inf)
  #   }
  #   val
  # }
  
  # do the lower bound
  h_low <- lp(x_low * Gam_k + theta_curr, dat)
  i <- 0
  n_expand <- 0
  while (i < max_steps_out & isTRUE(h_low > h)) {
    x_low <- x_low - width
    h_low <- lp(x_low * Gam_k + theta_curr, dat)
    i <- i + 1
  }
  int[1] <- x_low
  n_expand <- n_expand + i
  
  # now expand upper bound
  h_high <- lp(x_high * Gam_k + theta_curr, dat)
  i <- 0
  while(i < max_steps_out & isTRUE(h_high > h)) {
    x_high <- x_high + width
    h_high <- lp(x_high * Gam_k + theta_curr, dat)
    i <- i + 1
  }
  int[2] <- x_high
  n_expand <- n_expand + i
  
  if(n_expand == 0){ n_expand <- 1 }
  
  return(list(int = int, n_expand = n_expand))
  
}



#' Draw a proposal from a slice interval via shrinkage (contraction)
#'
#' Samples uniformly from the current interval \code{A_curr} along factor direction
#' \code{k} and contracts the interval toward the current point whenever the proposal
#' falls below the slice height \code{h}, following the shrinkage procedure of Neal (2003).
#'
#' @param k Integer index of the factor direction to sample along.
#' @param A_curr Numeric vector of length 2 giving the current lower and upper bounds
#'   of the slice interval on the standardised factor scale.
#' @param theta_curr Numeric vector of current parameter values.
#' @param pars_afs List of factor slice sampler parameters. Must contain
#'   \code{Gamma}, a matrix whose columns are the factor directions.
#' @param h Numeric scalar. Log-probability threshold defining the current slice height.
#' @param lp Function to evaluate the log-joint probability. Must accept a parameter
#'   vector as its first argument and the data list as its second.
#' @param dat Data list passed to \code{lp}.
#' @param time_out Maximum number of proposals before an error is thrown. Defaults
#'   to \code{100000}.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{theta}}{Numeric vector of the accepted proposal in parameter space.}
#'     \item{\code{contraction_count}}{Integer count of how many times the interval
#'       was contracted before a valid proposal was found.}
#'   }
#'
slice_proposals <- function(k, A_curr, theta_curr, pars_afs, h, lp, dat, time_out = 100000){
  
  h_prop <- -Inf
  prop_counter <- 0
  contraction_counter <- 0
  while(h_prop < h){
    q_prop <- runif(1, min = A_curr[1], max = A_curr[2])
    h_prop <- lp(q_prop * pars_afs$Gamma[, k] + theta_curr, dat)
    if(h_prop < h){
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
    stop(paste0("Failed to find suitable proposal for factor ", k))
  }
  theta_prop <- q_prop * pars_afs$Gamma[, k] + theta_curr
  
  return(
    list(theta = theta_prop, contraction_count = contraction_counter)
  )
  
}


#' Run the Factor Slice Sampler (Algorithm 2 of Tibbits et al. 2013)
#'
#' At each iteration a random slice height is drawn beneath the log-probability of
#' the current point, then each factor direction is updated in turn via
#' \code{\link{step_out}} (to bracket the slice) and \code{\link{slice_proposals}}
#' (to draw from within the bracket). The factor directions are taken from
#' \code{pars_afs$Gamma}, so the sampler reduces to univariate slice sampling along
#' the coordinate axes when \code{Gamma} is the identity matrix.
#'
#' @param theta_init Named numeric vector of initial parameter values.
#' @param dat Data list passed to \code{lp}.
#' @param lp Function to evaluate the log-joint probability. Must accept a parameter
#'   vector as its first argument and the data list as its second.
#' @param pars_afs List of factor slice sampler parameters. Must contain:
#'   \describe{
#'     \item{\code{Gamma}}{Matrix whose columns are the factor directions.}
#'     \item{\code{width}}{Numeric vector of interval widths, one per factor.}
#'   }
#' @param iter Number of MCMC iterations to run. Defaults to \code{1000}.
#' @param control Named list of algorithm controls. Recognised keys:
#'   \describe{
#'     \item{\code{time_out}}{Max proposals per contraction step (default \code{10000}).}
#'     \item{\code{max_steps_out}}{Max expansion steps per direction (default \code{10000}).}
#'     \item{\code{track_interval_steps}}{Logical; whether to record stepping-out and
#'       contraction counts (default \code{TRUE}).}
#'   }
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{theta}}{Numeric matrix with \code{iter} rows and one column per
#'       parameter; row names omitted, column names inherited from \code{theta_init}.}
#'     \item{\code{out_steps}}{Integer matrix (\code{iter} x \code{p}) of stepping-out
#'       counts per iteration and factor direction; \code{NULL} if
#'       \code{track_interval_steps} is \code{FALSE}.}
#'     \item{\code{in_steps}}{Integer matrix (\code{iter} x \code{p}) of contraction
#'       counts per iteration and factor direction; \code{NULL} if
#'       \code{track_interval_steps} is \code{FALSE}.}
#'   }
#'
factor_slice_sampler <- function(theta_init, dat, lp, pars_afs, iter = 1000, control = list()){
  
  # set controls for other funs if not set
  alg_control <- list(
    time_out = 10000,
    max_steps_out = 10000,
    track_interval_steps = TRUE
  )
  alg_control <- modifyList(alg_control, control)
  
  # initialize trackers
  p <- length(theta_init)
  theta <- matrix(nrow = iter + 1, ncol = p)
  colnames(theta) <- names(theta_init)
  theta[1, ] <- theta_init
  
  if(alg_control$track_interval_steps){
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
      samp_i <- slice_proposals(j, A_j, theta_prop, pars_afs, h_i, lp, dat, time_out = alg_control$time_out)
      
      theta_prop <- samp_i$theta
      
      if(alg_control$track_interval_steps){
        out_steps[i, j] <- out$n_expand
        in_steps[i, j] <- samp_i$contraction_count
      }
      
    }
    theta[i + 1, ] <- theta_prop
  }
  ret <- list(theta = theta[-1, ])
  if(alg_control$track_interval_steps){
    ret$out_steps <- out_steps
    ret$in_steps <- in_steps
  }
  return(ret)
}



#' Tune factor slice sampler interval widths
#'
#' Runs the factor slice sampler and iteratively adjusts the interval width for
#' each factor direction so that the ratio of stepping-out steps to total steps
#' (stepping-out + contraction) converges toward \code{ratio_target}. Widths are
#' scaled up when the ratio exceeds the target and scaled down when it falls below.
#' Tuning is declared converged for a given direction once the multiplicative
#' adjustment falls within 10\% of 1.
#'
#' @param theta_init Named numeric vector of initial parameter values.
#' @param dat Data list passed to \code{lp}.
#' @param lp Function to evaluate the log-joint probability. Must accept a parameter
#'   vector as its first argument and the data list as its second.
#' @param pars_afs List of factor slice sampler parameters, or \code{NULL} to
#'   initialise with identity factor directions and small random widths. Must contain
#'   \code{Gamma} (factor direction matrix) and \code{width} (per-factor widths)
#'   when not \code{NULL}.
#' @param tune_w_after Number of sampler steps between width adjustments. Doubles
#'   automatically once all step-out rates drop below 10 per adjustment window.
#'   Defaults to \code{1}.
#' @param ratio_target Target ratio of stepping-out steps to total steps (stepping-out
#'   + contraction). Defaults to \code{0.8}.
#' @param stop_after Maximum total sampler steps before halting, regardless of
#'   convergence. Defaults to \code{2^10 + 1}.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{pars_afs}}{Updated parameter list with tuned \code{width} values.}
#'     \item{\code{theta}}{Named numeric vector of the last sampled parameter values,
#'       suitable for use as \code{theta_init} in a subsequent run.}
#'   }
#'
tune_widths <- function(theta_init, dat, lp, pars_afs = NULL, init_width = 0.1, tune_w_after = 1, ratio_target = 0.8, stop_after = 2^10 + 1){
  
  p <- length(theta_init)
  
  # ---- Initialize sampler and conduct tuning ----
  if(is.null(pars_afs)){
    pars_afs <- list(
      Gamma = diag(nrow = p),
      width = runif(p, max = init_width),
      theta_cov = diag(nrow = p)
    )
  }
  theta <- theta_init
  
  # set counters
  iter <- 1
  counter_w <- 1
  # counter_cov <- 1
  convergence <- rep(1, p)
  names(convergence) <- paste("theta", 1:p, sep = "_")
  steps_out <- rep(0, p)
  contractions <- rep(0, p)
  while(iter <= stop_after & any(convergence == 1)){
    
    samp_i <- factor_slice_sampler(theta, dat, lp, pars_afs, iter = 1)
    
    # add to running totals
    steps_out <- samp_i$out_steps[1, ] + steps_out
    contractions <- samp_i$in_steps[1, ] + contractions
    
    theta <- samp_i$theta
    
    # xbar <- xbar + samp_i$theta
    # # note that theta is a row vector here
    # sample_cov <- sample_cov + samp_i$theta %*% t(samp_i$theta)
    
    ## ---- Tuning the widths ----
    if(counter_w == tune_w_after){
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
      # this is just a "good enough" sort of criterion
      if(all((steps_out / tune_w_after) < 10) & all(contractions < steps_out)){
        convergence <- rep(0, p)
      } else if(all((steps_out / tune_w_after) < 10)) {
        tune_w_after <- 2 * tune_w_after
      }
      # end if
    } # end if
    counter_w <- counter_w + 1
    iter <- iter + 1
    
  } # end while
  
  return(list(
    pars_afs = pars_afs,
    theta = theta
  ))
}




#' Estimate factor directions from a burn-in run
#'
#' Runs the factor slice sampler for \code{burnin} iterations, computes the
#' empirical covariance matrix of the resulting samples, and sets the factor
#' directions (\code{pars_afs$Gamma}) to the eigenvectors of that covariance
#' matrix. This aligns the slice directions with the principal axes of the
#' posterior, which is the core adaptation step of the Factor Slice Sampler
#' (Tibbits et al. 2013).
#'
#' @param theta_init Named numeric vector of initial parameter values.
#' @param dat Data list passed to \code{lp}.
#' @param lp Function to evaluate the log-joint probability. Must accept a parameter
#'   vector as its first argument and the data list as its second.
#' @param pars_afs List of current factor slice sampler parameters. Must contain
#'   \code{Gamma} and \code{width}. Both are passed to \code{\link{factor_slice_sampler}};
#'   only \code{Gamma} and \code{theta_cov} are updated in the returned object.
#' @param burnin Number of iterations to run for covariance estimation. Defaults
#'   to \code{500}.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{pars_afs}}{Updated parameter list with \code{Gamma} set to the
#'       eigenvectors of the sample covariance and \code{theta_cov} set to the
#'       sample covariance matrix itself.}
#'     \item{\code{theta}}{Named numeric vector of the last sampled parameter values
#'       from the burn-in run, suitable for use as \code{theta_init} in a subsequent
#'       run.}
#'   }
#'
estimate_factors <- function(theta_init, dat, lp, pars_afs, burnin = 500){
  
  theta_burn <- factor_slice_sampler(theta_init, dat, lp, pars_afs, iter = burnin)$theta
  # estimate the covariance
  V <- cov(theta_burn)
  
  # estimate factors
  pars_afs$Gamma <- eigen(V)$vectors
  pars_afs$theta_cov <- V
  
  return(list(
    pars_afs = pars_afs, 
    theta = theta_burn[burnin, ]
  ))
  
}


#' Default control parameters for \code{\link{auto_tune_sampler}}
#'
#' Returns a named list of default values for the tuning pipeline. Pass the
#' output (or a modified copy) as the \code{control} argument of
#' \code{\link{auto_tune_sampler}}.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{\code{pars_afs}}{\code{NULL}; triggers initialisation from the
#'       identity matrix and small random widths inside \code{\link{tune_widths}}.}
#'     \item{\code{tune_w_after}}{Integer. Number of sampler steps between each
#'       width adjustment (default \code{1}).}
#'     \item{\code{ratio_target}}{Numeric. Target ratio of stepping-out to total
#'       steps (default \code{0.8}).}
#'     \item{\code{stop_after}}{Integer. Maximum width-tuning steps per round
#'       (default \code{2^10 + 1}).}
#'     \item{\code{burnin}}{Integer. Burn-in iterations used by
#'       \code{\link{estimate_factors}} to estimate the posterior covariance
#'       (default \code{500}).}
#'   }
#'
control_tuning <- function(){
  list(
    pars_afs = NULL,
    init_width = 0.1,
    tune_w_after = 1,
    ratio_target = 0.8,
    stop_after = 2^10 + 1,
    burnin = 500
  )
}


#' Automatically tune the Factor Slice Sampler
#'
#' Runs a three-stage auto-tuning pipeline for the Factor Slice Sampler:
#' \enumerate{
#'   \item \strong{Width tuning (round 1):} \code{\link{tune_widths}} is called
#'     with coordinate-aligned factor directions to find good initial interval widths.
#'   \item \strong{Factor estimation:} \code{\link{estimate_factors}} is called to
#'     rotate the slice directions to align with the posterior covariance structure.
#'   \item \strong{Width re-tuning (round 2):} \code{\link{tune_widths}} is called
#'     again along the estimated factor directions to refine the widths.
#' }
#' The returned \code{pars_afs} object can be passed directly to
#' \code{\link{factor_slice_sampler}} for production sampling.
#'
#' @param theta_init Named numeric vector of initial parameter values.
#' @param dat Data list passed to \code{lp}.
#' @param lp Function to evaluate the log-joint probability. Must accept a parameter
#'   vector as its first argument and the data list as its second.
#' @param control Named list of tuning controls. Recognised keys (all optional):
#'   \describe{
#'     \item{\code{pars_afs}}{Initial \code{pars_afs} list; \code{NULL} to start
#'       from scratch (default \code{NULL}).}
#'     \item{\code{tune_w_after}}{Passed to \code{\link{tune_widths}} (default \code{1}).}
#'     \item{\code{ratio_target}}{Target step-out ratio; passed to
#'       \code{\link{tune_widths}} (default \code{0.8}).}
#'     \item{\code{stop_after}}{Maximum width-tuning steps; passed to
#'       \code{\link{tune_widths}} (default \code{2^10 + 1}).}
#'     \item{\code{burnin}}{Burn-in iterations for factor estimation; passed to
#'       \code{\link{estimate_factors}} (default \code{500}).}
#'   }
#'
#' @return A list with two elements (the output of the final \code{\link{tune_widths}}
#'   call):
#'   \describe{
#'     \item{\code{pars_afs}}{Fully tuned parameter list with optimised \code{Gamma}
#'       (factor directions), \code{width} (per-factor interval widths), and
#'       \code{theta_cov} (estimated posterior covariance).}
#'     \item{\code{theta}}{Named numeric vector of the last sampled parameter values,
#'       ready to be used as the starting point for production sampling.}
#'   }
#'
auto_tune_sampler <- function(theta_init, dat, lp, control = list()){
  
  alg_control <- modifyList(control_tuning(), control)
  
  all_pars <- c(
    list(
      theta_init = theta_init,
      dat = dat,
      lp = lp
    ),
    alg_control
  )
  
  round1 <- with(all_pars, {
      tune_widths(
        theta_init,
        dat,
        lp,
        pars_afs,
        init_width,
        tune_w_after,
        ratio_target,
        stop_after
      )
  })
  
  # now estimate factors
  names(round1)[which(names(round1) == "theta")] <- "theta_init"
  all_pars <- modifyList(all_pars, round1)
  
  round2 <- with(all_pars,{
    estimate_factors(
      theta_init,
      dat,
      lp,
      pars_afs,
      burnin
    )
  })
  
  # update the list
  names(round2)[which(names(round2) == "theta")] <- "theta_init"
  all_pars <- modifyList(all_pars, round2)
  
  # retune the step sizes
  round3 <- with(all_pars, {
    tune_widths(
      theta_init,
      dat,
      lp,
      pars_afs,
      tune_w_after,
      ratio_target,
      stop_after
    )
  })
  
  return(round3)
  
}


#' Default control parameters for \code{\link{fss_in_Gibbs_sample}}
#'
#' Returns a named list of default values for the factor slice sampler when used
#' inside a Gibbs loop. Pass the output (or a modified copy) as the
#' \code{control_sampler} argument of \code{\link{fss_in_Gibbs_sample}}.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{\code{iter}}{Integer. Total number of Gibbs iterations to run
#'       (default \code{1000}).}
#'     \item{\code{time_out}}{Integer. Maximum contraction attempts per factor
#'       direction before \code{\link{slice_proposals}} throws an error
#'       (default \code{10000}).}
#'     \item{\code{max_steps_out}}{Integer. Maximum expansion steps per direction
#'       in \code{\link{step_out}} (default \code{10000}).}
#'     \item{\code{track_interval_steps}}{Logical. Whether to record stepping-out
#'       and contraction counts during sampling (default \code{TRUE}).}
#'   }
#'
control_fssig <- function(){
  return(
    list(
      iter = 1000,
      time_out = 10000,
      max_steps_out = 10000,
      track_interval_steps = TRUE
    )
  )
}


#' Factor Slice Sampler within a Gibbs data-augmentation loop
#'
#' Alternates between two steps at each iteration:
#' \enumerate{
#'   \item \strong{Gibbs step:} \code{fill_rng} is called to draw the missing
#'     observations from their full conditional distribution given the current
#'     parameter values \eqn{\boldsymbol\theta}.
#'   \item \strong{Factor slice step:} \code{\link{factor_slice_sampler}} is called
#'     for a single iteration to draw a new \eqn{\boldsymbol\theta} given the
#'     just-imputed complete data.
#' }
#' This function assumes that the factor slice sampler has already been tuned
#' (via \code{\link{auto_tune_sampler}}) and that a valid \code{pars_afs} object
#' is supplied.
#'
#' @param theta_init Named numeric vector of initial parameter values.
#' @param dat Data list. Must contain a \code{y} element that is a data frame with
#'   columns \code{yt} and \code{ytm1}, with \code{NA} in \code{yt} for any missing
#'   time steps.
#' @param lp Function to evaluate the log-joint probability. Must accept a parameter
#'   vector as its first argument and the data list as its second.
#' @param fill_rng Function that imputes missing observations. Must accept the current
#'   parameter vector as its first argument and the data list as its second, and return
#'   an updated version of \code{dat$y} with all \code{NA}s in \code{yt} filled in.
#' @param pars_afs Tuned factor slice sampler parameter list, as returned by
#'   \code{\link{auto_tune_sampler}}. Must contain \code{Gamma} (factor direction
#'   matrix) and \code{width} (per-factor interval widths).
#' @param control_sampler Named list of sampler controls. Unrecognised keys are
#'   ignored; missing keys fall back to \code{\link{control_fssig}} defaults.
#'   Recognised keys:
#'   \describe{
#'     \item{\code{iter}}{Total number of Gibbs iterations (default \code{1000}).}
#'     \item{\code{time_out}}{Max contraction attempts per factor direction
#'       (default \code{10000}).}
#'     \item{\code{max_steps_out}}{Max expansion steps per direction
#'       (default \code{10000}).}
#'   }
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{theta}}{Numeric matrix with \code{iter} rows and one column per
#'       parameter; column names inherited from \code{theta_init}.}
#'     \item{\code{y}}{Numeric matrix with \code{iter} rows and one column per time
#'       step, storing the full (imputed) \code{yt} vector at each iteration.
#'       Columns corresponding to observed time steps will be constant.}
#'     \item{\code{miss_ids}}{Integer vector of row indices in \code{dat$y} where
#'       \code{yt} was \code{NA}, identifying which columns of \code{y} are imputed.}
#'   }
#'
fss_in_Gibbs_sample <- function(theta_init, dat, lp, fill_rng, pars_afs, control_sampler = list()){
  
  cntrl_fss <- modifyList(control_fssig(), control_sampler)
  
  iter <- cntrl_fss$iter
  p    <- length(theta_init)
  
  # storage for parameter samples
  theta <- matrix(nrow = iter, ncol = p)
  colnames(theta) <- names(theta_init)
  
  # track the full yt column so callers can inspect imputed values
  y_samps  <- matrix(nrow = iter, ncol = nrow(dat$y))
  miss_ids <- which(is.na(dat$y[, "yt"]))
  
  # inner FSS control: one slice step per Gibbs cycle, no step tracking
  fss_ctrl <- list(
    time_out = cntrl_fss$time_out,
    max_steps_out = cntrl_fss$max_steps_out,
    track_interval_steps = FALSE
  )
  
  # recycled variables
  theta_curr <- theta_init
  dat_s <- dat
  
  for (i in 1:iter) {
    
    # Gibbs step: draw missing y given current theta 
    dat_s$y <- fill_rng(theta_curr, dat)
    
    # Factor slice step: draw theta given filled-in data
    fss_out <- factor_slice_sampler(
      theta_curr, 
      dat_s, 
      lp, 
      pars_afs,
      iter = 1,
      control = fss_ctrl
    )
    # update theta
    theta_curr <- fss_out$theta
    
    # Store
    theta[i, ]   <- theta_curr
    y_samps[i, ] <- dat_s$y[, "yt"]
    
  }
  
  return(list(
    theta = theta,
    y = y_samps,
    miss_ids = miss_ids
  ))
  
  
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
#' Estimates the posterior distribution of the Ricker model parameters \eqn{r}
#' and \eqn{\alpha} (and dispersion \eqn{\psi} for \code{fam = "neg_binom"}) using
#' a data-augmentation MCMC scheme. Starting values are obtained via empirical Bayes
#' (complete-case MLE), the Factor Slice Sampler is auto-tuned on the complete cases,
#' and then \code{chains} independent chains are run via
#' \code{\link{fss_in_Gibbs_sample}}, which alternates between Gibbs imputation of
#' missing counts and factor slice updates of the parameters.
#'
#' @param y Observed count data, supplied in one of two formats depending on
#'   \code{off_patch}:
#'   \describe{
#'     \item{Vector (\code{off_patch = FALSE})}{A numeric vector of counts in
#'       chronological order. \code{NA} indicates a missing observation. Any
#'       leading \code{NA}s are stripped with a warning because the first
#'       observation is required as a known initial condition.}
#'     \item{Data frame (\code{off_patch = TRUE})}{A data frame already in
#'       lagged-pair format with columns \code{yt} (count at time \eqn{t}) and
#'       \code{ytm1} (count at time \eqn{t-1}). \code{NA} in \code{yt} indicates
#'       a missing observation. Use this format when observations come from
#'       multiple patches or when the lagged pairs have been pre-constructed
#'       externally.}
#'   }
#' @param fam Character string specifying the observation model. One of
#'   \code{"poisson"} (default) or \code{"neg_binom"}.
#' @param chains Positive integer. Number of independent MCMC chains to run. Each
#'   produces \code{samples} draws from the posterior. Default \code{3}.
#' @param samples Positive integer. Number of posterior samples to retain per
#'   chain. Default \code{1000}.
#' @param priors_list Named list of prior hyperparameters. All priors are normal.
#'   Recognised keys (defaults in parentheses):
#'   \describe{
#'     \item{\code{m_r}, \code{sd_r}}{Mean and SD for the prior on \eqn{r}
#'       (\code{0}, \code{2.5}).}
#'     \item{\code{m_lalpha}, \code{sd_lalpha}}{Mean and SD for the prior on
#'       \eqn{\log\alpha} (\code{-3}, \code{1}).}
#'     \item{\code{m_lpsi}, \code{sd_lpsi}}{Mean and SD for the prior on
#'       \eqn{\log\psi}; only used when \code{fam = "neg_binom"}
#'       (\code{0}, \code{100}).}
#'   }
#' @param nthin Positive integer. Thinning interval; only every \code{nthin}-th
#'   sample is retained. Default \code{1} (no thinning).
#' @param return_y Logical. If \code{TRUE} and there are missing observations,
#'   the posterior samples of each imputed count are appended as extra columns
#'   (named \code{y1}, \ldots, \code{yn}) to all summary outputs. Columns
#'   corresponding to observed time steps will be constant. Default \code{FALSE}.
#' @param off_patch Logical. If \code{TRUE}, \code{y} is expected to be a data
#'   frame with columns \code{yt} and \code{ytm1} (already in lagged-pair format,
#'   e.g. for a multi-patch study). If \code{FALSE} (default), \code{y} is a plain
#'   vector and the lagged-pair data frame is constructed internally.
#' @param control_sampler Named list of controls passed to
#'   \code{\link{fss_in_Gibbs_sample}} via \code{\link{control_fssig}}. Use this
#'   to override the number of sampling iterations, time-out limits, etc.
#' @param control_tuner Named list of controls passed to
#'   \code{\link{auto_tune_sampler}} via \code{\link{control_tuning}}. Use this
#'   to override the width-tuning schedule, factor-estimation burn-in, etc.
#'
#' @return On success, a named list with five elements:
#'   \describe{
#'     \item{\code{estim}}{Named numeric vector of posterior means. Parameters are
#'       returned on the natural scale (\eqn{\alpha}, not \eqn{\log\alpha};
#'       \eqn{\psi}, not \eqn{\log\psi}). If \code{return_y = TRUE} and data are
#'       missing, imputed count means are appended as \code{y1}, \ldots, \code{yn}.}
#'     \item{\code{rhat}}{Named numeric vector of Gelman-Rubin \eqn{\hat{R}}
#'       convergence diagnostics (one per parameter), computed via
#'       \code{posterior::rhat}.}
#'     \item{\code{se}}{Named numeric vector of posterior standard deviations.}
#'     \item{\code{lower}}{Named numeric vector of 2.5\% posterior quantiles.}
#'     \item{\code{upper}}{Named numeric vector of 97.5\% posterior quantiles.}
#'   }
#'   On failure (population extinction, insufficient non-missing pairs, sampler
#'   getting stuck, etc.) a two-element list \code{list(NA, reason = "...")} is
#'   returned with a descriptive warning.
#'
#' @seealso \code{\link{auto_tune_sampler}}, \code{\link{fss_in_Gibbs_sample}},
#'   \code{\link{control_tuning}}, \code{\link{control_fssig}}
#'
fit_ricker_DA <- function(
    y, fam = "poisson", 
    chains = 3, 
    priors_list = list(
      m_r = 0,
      sd_r = 2.5,
      m_lalpha = -3,
      sd_lalpha = 1,
      m_lpsi = 0,
      sd_lpsi = 100
    ),
    nthin = 1,
    return_y = FALSE,
    off_patch = FALSE,
    control_sampler = list(),
    control_tuner = list()
){

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
  writeLines("Finding intial parameter values using complete cases...\n")
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
  names(theta_init)[which(names(theta_init) == "alpha")] <- "lalpha"
  
  
  # some tuning if fam = neg_binom
  if(fam == "neg_binom"){
    if(is.infinite(theta_init["psi"])){
      theta_init["psi"] <- log(100)
    } else {
      theta_init["psi"] <- log(theta_init["psi"])
    }
    names(theta_init)[which(names(theta_init) == "psi")] <- "lpsi"
  }
  
  
  
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
    # defend against infinite or NaN
    if(!is.finite(lprob) | is.nan(lprob) | is.na(lprob)){
      return(-Inf)
    }
    return(unname(lprob))
  }
  
  # function to "fill in" the missing data based on current values of params
  fill_rng <- function(theta, datlist){
    y <- datlist$y
    for(t in 1:nrow(y)){
      if(is.na(y[t, "yt"])){
        mu_t <- ricker_step(
          theta = c(theta["r"], exp(theta["lalpha"])), 
          Nt = y[t, "ytm1"]
        )
        # fallback
        if(!is.finite(mu_t)){
          warning("computed non finite mu_t. Falling back to y_tm1")
          mu_t <- y[t, "ytm1"]
        }
        if(datlist$fam == "poisson"){
          y[t, "yt"] <- zt_poisson_rng(1, mu_t)
          if(is.na(y[t + 1, "ytm1"])){
            y[t + 1, "ytm1"] <- y[t, "yt"]
          }
        }
        if(datlist$fam == "neg_binom"){
          y[t, "yt"] <- zt_neg_binom_rng(1, size = exp(theta["lpsi"]), mu = mu_t)
          if(is.na(y[t + 1, "ytm1"])){
            y[t + 1, "ytm1"] <- y[t, "yt"]
          }
        }
      }
    }
    return(y)
  }
  
  # compile data
  datlist <- c(
    list(
      y = dat,
      fam = fam
    ),
    priors_list
  )
  
  datlist_cc <- datlist
  datlist_cc$y <- dat[complete.cases(dat), ]
  
  # ---- Conduct tuning ----
  writeLines("Tuning the sampler...\n")
  
  sampler_tune <- auto_tune_sampler(
    theta_init, 
    dat = datlist_cc,
    lp = compute_lp,
    control = control_tuner
  )
  
  # ---- Main sampling ----
  writeLines("Sampling from posterior...\n")
  post_samps <- replicate(
    chains,
    fss_in_Gibbs_sample(
      theta_init = sampler_tune$theta,
      dat = datlist,
      lp = compute_lp,
      fill_rng = fill_rng,
      pars_afs = sampler_tune$pars_afs,
      control_sampler = control_sampler
    ),
    simplify = FALSE
  )
  
  # compute rhat statistic for samples
  samps_by_chains <- lapply(
    names(theta_init),
    FUN = function(n, post_samps){
      Reduce(
        cbind,
        lapply(post_samps, \(x) x$theta[, n] )
      )
    },
    post_samps = post_samps
  )
  
  rhats <- sapply(samps_by_chains, posterior::rhat)
  
  ses <- apply(
    Reduce(cbind, samps_by_chains),
    2, 
    sd
  )
  names(ses) <- as.vector(sapply(
    names(theta_init),
    \(n) paste0(n, 1:chains)
  ))
  if(any(ses == 0)){
    stuck <- names(ses)[which(ses == 0)]
    mess <- paste(paste(
      "MCMC sampler got stuck. Chain", stuck
    ), collapse = "\n")
    warning(mess)
    return(list(
      NA,
      reason = "Sampler got stuck"
    ))
  }
  
  # combine posterior samps of r and alpha
  theta_samps <- Reduce(
    rbind,
    lapply(post_samps, \(x){ x$theta })
  )
  theta_samps[,"lalpha"] <- exp(theta_samps[,"lalpha"])
  colnames(theta_samps)[1:2] <- c("r", "alpha")
  if(fam == "neg_binom"){
    theta_samps[, "lpsi"] <- exp(theta_samps[, "lpsi"])
    colnames(theta_samps)[3] <- "psi"
  }
  
  # if we want to return posterior estims for missing obs
  prop_miss <- mean(is.na(datlist$y[, "yt"]))
  if(isTRUE(return_y) & prop_miss > 0){
    y_samps <- Reduce(
      rbind,
      lapply(post_samps, function(x){
        x$y
      })
    )
    names_all <- c(colnames(theta_samps), paste0("y", 1:n))
    theta_samps <- cbind(
      theta_samps, y_samps
    )
    colnames(theta_samps) <- names_all
  }
  
  # return summaries
  return(
    list(
      estim = apply(theta_samps, 2, mean),
      rhat = rhats,
      se = apply(theta_samps, 2, sd),
      lower = apply(theta_samps, 2, quantile, probs = 0.025),
      upper = apply(theta_samps, 2, quantile, probs = 0.975)
    )
  )
  
}

