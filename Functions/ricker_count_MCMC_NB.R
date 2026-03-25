###############################################
# MCMC samplers for Bayesian Ricker model
# Supports:
#   - Single or multiple independent time series
#   - Poisson or Negative Binomial error
#   - Missing data via Data Augmentation
###############################################


#' Compute the mode of an empirical distribution
#'
#' @param x Vector of samples from the target distribution.
#' @return The mode of the empirical distribution.
#'
posterior_mode <- function(x){
  d <- density(x)
  return(d$x[which.max(d$y)])
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Build a model matrix X for a single complete y vector
.make_X <- function(y){
  n <- length(y)
  cbind(rep(1, n), y)
}

#' Initialize theta from a list of y vectors
#'
#' Fits \code{fit_ricker_cc} to each series, averages estimates, and
#' returns a named starting vector. For neg_binom adds lpsi.
#'
.init_theta <- function(y_list, fam = "poisson"){

  fits <- lapply(y_list, function(y){
    fit_ricker_cc(y, fam = fam)
  })

  # drop failed fits
  ok <- sapply(fits, function(f) !is.null(f$estim) && length(f$estim) > 1 && !is.na(f$estim[1]))
  fits <- fits[ok]

  if(length(fits) == 0)
    stop("fit_ricker_cc failed on all series during initialization.")

  r_vals     <- sapply(fits, function(f) f$estim[1])
  alpha_vals <- sapply(fits, function(f) f$estim[2])

  r_mean     <- mean(r_vals,     na.rm = TRUE)
  alpha_mean <- mean(alpha_vals, na.rm = TRUE)

  # alpha from fit_ricker_cc is already positive (negated internally);
  # convert to log scale for the sampler
  if(alpha_mean <= 0){
    # fall back: use the upper CI average
    alpha_upper <- sapply(fits, function(f) f$upper[2])
    alpha_mean  <- mean(alpha_upper, na.rm = TRUE)
    if(alpha_mean <= 0) alpha_mean <- 0.001
  }
  lalpha <- log(alpha_mean)

  theta_init <- c(r = r_mean, lalpha = lalpha)

  if(fam == "neg_binom"){
    # extract dispersion from MASS::glm.nb fits where available
    psi_vals <- sapply(fits, function(f){
      if(!is.null(f$theta)) f$theta else NA_real_
    })
    psi_mean <- mean(psi_vals, na.rm = TRUE)
    if(is.na(psi_mean) || psi_mean <= 0) psi_mean <- 1.0
    theta_init <- c(theta_init, lpsi = log(psi_mean))
  }

  return(theta_init)
}

#' Compute the Hessian-based proposal covariance for a list of series
#'
#' Uses the joint NLL across all series evaluated at \code{theta_init}.
#'
.proposal_sigma <- function(theta_init, y_list, fam = "poisson"){

  # build complete-case X matrices
  X_list <- lapply(y_list, .make_X)

  # NLL wrapper that works on natural-scale psi for neg_binom
  nll_joint <- function(x){
    theta_nll <- x
    if(fam == "neg_binom"){
      # x[3] is lpsi; convert to psi for the likelihood
      theta_nll[3] <- exp(x[3])
    }
    ricker_count_neg_ll_multi(theta_nll, y_list, X_list, fam)
  }

  hess <- tryCatch(
    optim(theta_init, nll_joint, hessian = TRUE)$hessian,
    error = function(e){
      warning("Hessian computation failed, using diagonal fallback: ", e$message)
      diag(0.1, length(theta_init))
    }
  )

  # regularize if near-singular
  eigs <- eigen(hess, only.values = TRUE)$values
  if(any(eigs < 1e-8)){
    warning("Hessian near-singular, regularizing.")
    hess <- hess + diag(1e-6, nrow(hess))
  }

  Sigma <- tryCatch(solve(hess), error = function(e){
    warning("solve(hess) failed, using diag fallback.")
    diag(0.1, nrow(hess))
  })

  return(Sigma)
}


# ---------------------------------------------------------------------------
# Proposal distribution (shared by both samplers)
# ---------------------------------------------------------------------------

.q_rng <- function(theta, Sigma){
  prop <- mvtnorm::rmvnorm(1, theta, Sigma)[1, ]
  names(prop) <- names(theta)
  prop
}

.q_lpdf <- function(prop, theta, Sigma){
  mvtnorm::dmvnorm(prop, theta, Sigma, log = TRUE)
}


# ---------------------------------------------------------------------------
# MH block sampler (no missing data)
# ---------------------------------------------------------------------------

#' Metropolis-Hastings sampling of the posterior with block updates
#'
#' @param dat List with elements \code{y_list} (list of complete count vectors),
#'   plus prior hyperparameters \code{m_r}, \code{sd_r}, \code{m_lalpha},
#'   \code{sd_lalpha}, and (for neg_binom) \code{m_lpsi}, \code{sd_lpsi}.
#' @param lp Log-joint probability function: \code{lp(theta, dat)}.
#' @param burnin Number of warm-up samples.
#' @param iter Number of samples to keep.
#' @param nthin Thinning ratio.
#'
#' @return List with \code{theta} matrix and \code{acceptance} vector.
#'
MH_block_sample <- function(dat, lp, burnin, iter, nthin = 1){

  theta_init <- .init_theta(dat$y_list, fam = dat$fam)
  Sigma      <- .proposal_sigma(theta_init, dat$y_list, fam = dat$fam)

  S <- (burnin + iter) * nthin
  theta_samps <- matrix(nrow = S, ncol = length(theta_init))
  colnames(theta_samps) <- names(theta_init)
  theta_samps[1, ] <- theta_init

  lp_samps <- vector("double", S)
  lp_samps[1] <- lp(theta_init, dat)
  accept <- vector("double", S)

  for(s in 2:S){
    theta_s    <- theta_samps[s - 1, ]
    theta_prop <- .q_rng(theta_s, Sigma)
    lp_prop    <- lp(theta_prop, dat)
    lp_curr    <- lp_samps[s - 1]

    mh_ratio <- exp(
      (lp_prop + .q_lpdf(theta_s, theta_prop, Sigma)) -
      (lp_curr + .q_lpdf(theta_prop, theta_s, Sigma))
    )
    if(is.nan(mh_ratio)) mh_ratio <- 0

    A        <- min(1, mh_ratio)
    accept[s] <- rbinom(1, 1, prob = A)

    if(accept[s] == 1){
      theta_samps[s, ] <- theta_prop
      lp_samps[s]      <- lp_prop
    } else {
      theta_samps[s, ] <- theta_s
      lp_samps[s]      <- lp_curr
    }
  }

  samps2keep <- seq(1, S, by = nthin)
  theta_out  <- theta_samps[samps2keep, ][(burnin + 1):(burnin + iter), ]
  colnames(theta_out) <- names(theta_init)

  list(theta = theta_out, acceptance = accept)
}


# ---------------------------------------------------------------------------
# MH-within-Gibbs with Data Augmentation
# ---------------------------------------------------------------------------

#' Block Metropolis-Hastings within Gibbs sampling with Data Augmentation
#'
#' @param dat List with \code{y_list} (list of count vectors, NAs allowed),
#'   prior hyperparameters, and \code{fam}.
#' @param fill_rng Function to impute missing values: \code{fill_rng(theta, dat)},
#'   returns an updated \code{y_list} with NAs replaced.
#' @param lp Log-joint probability: \code{lp(theta, dat, y_list_full)}.
#' @param burnin Number of warm-up samples.
#' @param iter Number of samples to keep.
#' @param nthin Thinning ratio.
#'
#' @return List with \code{theta} matrix, \code{y_list} (list of matrices,
#'   one row per kept sample), and \code{acceptance} vector.
#'
MH_Gibbs_DA <- function(dat, fill_rng, lp, burnin, iter, nthin = 1){

  S <- (burnin + iter) * nthin

  # --- initialize ---
  theta_init    <- .init_theta(dat$y_list, fam = dat$fam)
  y_full_init   <- fill_rng(theta_init, dat)        # list of imputed y vectors
  Sigma         <- .proposal_sigma(theta_init, y_full_init, fam = dat$fam)

  # storage
  theta_samps <- matrix(nrow = S, ncol = length(theta_init))
  colnames(theta_samps) <- names(theta_init)
  theta_samps[1, ] <- theta_init

  # y stored as a list-of-matrices (one matrix per series, rows = iterations)
  n_series <- length(dat$y_list)
  y_samps  <- lapply(seq_len(n_series), function(i){
    matrix(nrow = S, ncol = length(dat$y_list[[i]]))
  })
  for(i in seq_len(n_series)) y_samps[[i]][1, ] <- y_full_init[[i]]

  lp_samps <- vector("double", S)
  lp_samps[1] <- lp(theta_init, dat, y_full_init)

  # if initial lp is -Inf, try random restarts
  iter_restart <- 1
  while(is.infinite(lp_samps[1]) && iter_restart < 100){
    theta_rt <- runif(length(theta_init), -1, 1)
    names(theta_rt) <- names(theta_init)
    y_try <- fill_rng(theta_rt, dat)
    lp_samps[1] <- lp(theta_rt, dat, y_try)
    if(!is.infinite(lp_samps[1])){
      theta_samps[1, ] <- theta_rt
      for(i in seq_len(n_series)) y_samps[[i]][1, ] <- y_try[[i]]
    }
    iter_restart <- iter_restart + 1
  }
  if(is.infinite(lp_samps[1]))
    stop("Failed to find suitable starting values.")

  accept  <- vector("double", S)
  accept2 <- vector("double", S)

  for(s in 2:S){

    # --- step 1: update theta (MH) ---
    theta_s    <- theta_samps[s - 1, ]
    theta_prop <- .q_rng(theta_s, Sigma)

    # reconstruct current y_list from storage
    y_curr <- lapply(seq_len(n_series), function(i) y_samps[[i]][s - 1, ])

    dat_s        <- dat
    dat_s$y_list <- y_curr

    lp_prop <- lp(theta_prop, dat_s)
    lp_curr <- lp_samps[s - 1]

    mh_ratio <- exp(
      (lp_prop + .q_lpdf(theta_s, theta_prop, Sigma)) -
      (lp_curr + .q_lpdf(theta_prop, theta_s, Sigma))
    )
    if(is.nan(mh_ratio)) mh_ratio <- 0

    A         <- min(1, mh_ratio)
    accept[s] <- rbinom(1, 1, prob = A)

    if(accept[s] == 1){
      theta_samps[s, ] <- theta_prop
      lp_samps[s]      <- lp_prop
    } else {
      theta_samps[s, ] <- theta_s
      lp_samps[s]      <- lp_curr
    }

    # --- step 2: update missing y values (MH on imputation) ---
    theta_cur    <- theta_samps[s, ]
    y_full_prop  <- fill_rng(theta_cur, dat)
    lp_prop2     <- lp(theta_cur, dat, y_full_prop)
    lp_curr2     <- lp_samps[s]

    mh_ratio2 <- exp(lp_prop2 - lp_curr2)
    if(is.nan(mh_ratio2)) mh_ratio2 <- 0

    A2         <- min(1, mh_ratio2)
    accept2[s] <- rbinom(1, 1, prob = A2)

    if(accept2[s] == 1){
      for(i in seq_len(n_series)) y_samps[[i]][s, ] <- y_full_prop[[i]]
      lp_samps[s] <- lp_prop2
    } else {
      for(i in seq_len(n_series)) y_samps[[i]][s, ] <- y_curr[[i]]
    }
  }

  # thin and discard burnin
  samps2keep    <- seq(1, S, by = nthin)
  keep_idx      <- samps2keep[(burnin + 1):(burnin + iter)]

  theta_out <- theta_samps[keep_idx, ]
  colnames(theta_out) <- names(theta_init)

  y_out <- lapply(y_samps, function(m) m[keep_idx, ])

  list(
    theta      = theta_out,
    y_list     = y_out,
    acceptance = accept
  )
}


# ---------------------------------------------------------------------------
# Main user-facing wrapper
# ---------------------------------------------------------------------------

#' Fit a Bayesian Ricker population growth model to one or more time series,
#' with optional missing data imputation and negative binomial errors.
#'
#' @param y Either a single numeric vector of counts, or a \strong{list} of
#'   such vectors. \code{NA} marks missing observations. Each series must
#'   begin with a non-NA value.
#' @param fam Character: \code{"poisson"} or \code{"neg_binom"}.
#' @param chains Number of parallel MCMC chains.
#' @param samples Number of posterior samples to keep per chain.
#' @param burnin Number of warm-up samples per chain.
#' @param priors_list Named list of prior hyperparameters:
#'   \describe{
#'     \item{m_r, sd_r}{Normal prior on r}
#'     \item{m_lalpha, sd_lalpha}{Normal prior on log(alpha)}
#'     \item{m_lpsi, sd_lpsi}{Normal prior on log(psi); only used for neg_binom}
#'   }
#' @param nthin Thinning interval.
#' @param return_y Logical. If \code{TRUE} and there is missingness, return
#'   posterior summaries for imputed observations.
#'
#' @return A list with elements \code{estim}, \code{rhat}, \code{se},
#'   \code{lower}, \code{upper}. If \code{return_y = TRUE} and data have
#'   missing values, imputed observation summaries are appended.
#'
fit_ricker_DA <- function(
    y,
    fam      = "poisson",
    chains   = 4,
    samples  = 1000,
    burnin   = 5000,
    priors_list = list(
      m_r       = 0,    sd_r       = 2.5,
      m_lalpha  = -3,   sd_lalpha  = 1,
      m_lpsi    = 0,    sd_lpsi    = 1     # used only for neg_binom
    ),
    nthin    = 5,
    return_y = FALSE
){
  require(parallel)

  # ------------------------------------------------------------------
  # 1. Coerce input to a list of y vectors
  # ------------------------------------------------------------------
  if(is.numeric(y)) y_list <- list(y)
  else if(is.list(y)) y_list <- y
  else stop("`y` must be a numeric vector or a list of numeric vectors.")

  n_series <- length(y_list)

  # ------------------------------------------------------------------
  # 2. Per-series validation
  # ------------------------------------------------------------------
  valid <- vector("list", n_series)
  for(i in seq_len(n_series)){
    yi <- y_list[[i]]

    if(any(yi == 0, na.rm = TRUE)){
      warning(sprintf("Series %d: population extinction (zeros found), returning NA.", i))
      return(list(NA, cause = "population extinction"))
    }
    if(any(is.nan(yi), na.rm = TRUE)){
      warning(sprintf("Series %d: NaN found, returning NA.", i))
      return(list(NA, reason = "NaN found"))
    }
    if(any(is.infinite(yi), na.rm = TRUE)){
      warning(sprintf("Series %d: infinite value found, returning NA.", i))
      return(list(NA, reason = "population explosion"))
    }

    # strip leading NAs
    if(is.na(yi[1])){
      warning(sprintf("Series %d: removing leading NAs.", i))
      start <- min(which(!is.na(yi)))
      yi <- yi[start:length(yi)]
    }

    # check enough complete pairs
    n_i   <- length(yi)
    dat_i <- data.frame(yt = yi[2:n_i], ytm1 = yi[1:(n_i - 1)])
    if(nrow(dat_i[complete.cases(dat_i), ]) < 5){
      warning(sprintf("Series %d: too few complete observation pairs.", i))
      return(list(NA, reason = "missingness limit"))
    }

    y_list[[i]] <- yi
  }

  # ------------------------------------------------------------------
  # 3. Build dat list (shared across chains)
  # ------------------------------------------------------------------
  prop_miss <- mean(sapply(y_list, function(yi) mean(is.na(yi))))

  dat <- c(
    list(y_list = y_list, fam = fam),
    priors_list
  )

  # ------------------------------------------------------------------
  # 4. Define log-posterior (joint across all series)
  # ------------------------------------------------------------------
  compute_lp <- function(theta, dat, y_list_full = NULL){

    fam_    <- dat$fam
    y_use   <- if(!is.null(y_list_full)) y_list_full else dat$y_list

    r      <- theta["r"]
    alpha  <- -exp(theta["lalpha"])   # enforce negativity
    theta_ll <- c(r, alpha)

    if(fam_ == "neg_binom"){
      psi      <- exp(theta["lpsi"])
      theta_ll <- c(theta_ll, psi)
    }

    # joint log-likelihood across all series
    X_list <- lapply(y_use, .make_X)
    ll     <- -ricker_count_neg_ll_multi(theta_ll, y_use, X_list, fam_)

    # priors
    lp <- ll +
      dnorm(r,              mean = dat$m_r,      sd = dat$sd_r,      log = TRUE) +
      dnorm(theta["lalpha"], mean = dat$m_lalpha, sd = dat$sd_lalpha, log = TRUE)

    if(fam_ == "neg_binom"){
      lp <- lp + dnorm(theta["lpsi"], mean = dat$m_lpsi, sd = dat$sd_lpsi, log = TRUE)
    }

    return(unname(lp))
  }

  # ------------------------------------------------------------------
  # 5. Define fill_rng (impute missing values in all series)
  # ------------------------------------------------------------------
  fill_rng <- function(theta, dat){

    fam_   <- dat$fam
    r      <- theta["r"]
    alpha  <- -exp(theta["lalpha"])

    if(fam_ == "neg_binom") psi <- exp(theta["lpsi"])

    lapply(dat$y_list, function(yi){
      n_i    <- length(yi)
      y_full <- numeric(n_i)
      y_full[1] <- yi[1]

      for(t in 2:n_i){
        if(is.na(yi[t])){
          mu_t <- ricker_step(c(r, alpha), y_full[t - 1])
          y_full[t] <- if(fam_ == "poisson"){
            rpois(1, mu_t)
          } else {
            rnbinom(1, size = psi, mu = mu_t)
          }
        } else {
          y_full[t] <- yi[t]
        }
      }
      # guard against zeros (avoid log(0))
      y_full[y_full == 0] <- 1L
      y_full
    })
  }

  # ------------------------------------------------------------------
  # 6. Run chains in parallel
  # ------------------------------------------------------------------
  force(list(compute_lp, fill_rng, dat))
  cl <- parallel::makeCluster(chains)

  parallel::clusterEvalQ(cl, {
    source(here::here("Functions/ricker_drop_function.R"))
    source(here::here("Functions/ricker_count_MCMC_NB.R"))
    source(here::here("Functions/ricker_count_likelihood_functions_NB.R"))
  })

  parallel::clusterExport(
    cl,
    varlist = ls(envir = environment()),
    envir   = environment()
  )

  if(prop_miss == 0){
    post_samps <- parallel::clusterCall(
      cl, MH_block_sample,
      dat    = dat,
      lp     = compute_lp,
      burnin = burnin,
      iter   = samples,
      nthin  = nthin
    )
  } else {
    post_samps <- parallel::clusterCall(
      cl, MH_Gibbs_DA,
      dat      = dat,
      lp       = compute_lp,
      fill_rng = fill_rng,
      burnin   = burnin,
      iter     = samples,
      nthin    = nthin
    )
  }

  parallel::stopCluster(cl)

  # ------------------------------------------------------------------
  # 7. Convergence check
  # ------------------------------------------------------------------
  post_r      <- Reduce(cbind, lapply(post_samps, function(x) x$theta[, "r"]))
  post_lalpha <- Reduce(cbind, lapply(post_samps, function(x) x$theta[, "lalpha"]))

  check_stuck <- function(mat, name){
    ses <- apply(mat, 2, sd)
    if(any(ses == 0)){
      stuck <- which(ses == 0)
      warning(sprintf("MCMC sampler stuck for param %s in chain(s) %s.",
                      name, paste(stuck, collapse = ", ")))
      return(TRUE)
    }
    FALSE
  }

  if(check_stuck(post_r, "r") || check_stuck(post_lalpha, "lalpha")){
    return(list(NA, reason = "Sampler got stuck"))
  }

  if(fam == "neg_binom"){
    post_lpsi <- Reduce(cbind, lapply(post_samps, function(x) x$theta[, "lpsi"]))
    if(check_stuck(post_lpsi, "lpsi"))
      return(list(NA, reason = "Sampler got stuck"))
  }

  # ------------------------------------------------------------------
  # 8. Combine chains and back-transform
  # ------------------------------------------------------------------
  theta_samps <- Reduce(rbind, lapply(post_samps, function(x) x$theta))

  # back-transform log-scale parameters
  theta_samps[, "lalpha"] <- exp(theta_samps[, "lalpha"])
  colnames(theta_samps)[colnames(theta_samps) == "lalpha"] <- "alpha"

  if(fam == "neg_binom"){
    theta_samps[, "lpsi"] <- exp(theta_samps[, "lpsi"])
    colnames(theta_samps)[colnames(theta_samps) == "lpsi"] <- "psi"
  }

  # ------------------------------------------------------------------
  # 9. Optionally append imputed y posteriors
  # ------------------------------------------------------------------
  if(isTRUE(return_y) && prop_miss > 0){
    for(ser_i in seq_len(n_series)){
      y_ser_samps <- Reduce(
        rbind,
        lapply(post_samps, function(x) x$y_list[[ser_i]])
      )
      colnames(y_ser_samps) <- paste0("y", ser_i, "_t", seq_len(ncol(y_ser_samps)))
      theta_samps <- cbind(theta_samps, y_ser_samps)
    }
  }

  # ------------------------------------------------------------------
  # 10. Return summaries
  # ------------------------------------------------------------------
  rhat_vals <- c(
    r     = posterior::rhat(post_r),
    alpha = posterior::rhat(post_lalpha)
  )
  if(fam == "neg_binom") rhat_vals["psi"] <- posterior::rhat(post_lpsi)

  list(
    estim = apply(theta_samps, 2, mean),
    rhat  = rhat_vals,
    se    = apply(theta_samps, 2, sd),
    lower = apply(theta_samps, 2, quantile, probs = 0.025),
    upper = apply(theta_samps, 2, quantile, probs = 0.975)
  )
}
