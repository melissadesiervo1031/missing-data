

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
fit_ricker_cc <- function(y, fam = "poisson", pro_conf="none"){
  
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
  
  n <- length(y)
  
  
  # compile into sliced dataframe
  dat <- data.frame(
    yt = y[2:n],
    ytm1 = y[1:(n - 1)]
  )
  
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
  
  # fit with poisson
  if(fam == "poisson"){
    fit <- glm(
      yt ~ ytm1, data = dat_cc, 
      family = poisson,
      offset = log(ytm1)
    )
  }
  
  # or fit with negbinom
  if(fam == "neg_binom"){
    fit <- MASS::glm.nb(
      yt ~ ytm1 + offset(log(ytm1)), 
      data = dat_cc
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
fit_ricker_drop <- function(y, fam = "poisson", pro_conf="none"){
  
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
  
  n <- length(y)
  
  # compile into sliced dataframe
  dat <- data.frame(
    yt = y[2:n],
    ytm1 = y[1:(n - 1)]
  )
  
  # fit with poisson
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





