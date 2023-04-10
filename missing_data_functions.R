#/////////////////////
# Function to introduce missingness into time series data
# Alice Stears
# updated 17 February 2023
#/////////////////////
# load packages -----------------------------------------------------------
library(tidyverse)
# function ----------------------------------------------------------------
### goals: 
## remove data at increasing levels of missingness...
makeMissing <- function(timeSeries, # a time series in vector format (a single vector, not a list)         
                        typeMissing, # a character string indicating how missingness is to be          
                        # generated ("random", "evenChunks", "randChunks", or "minMax")         
                        propMiss = NULL, # an optional vector argument that gives the proportion(s)          
                        # of the data that should be missing; if not specified, then increasingly          
                        # more missingness is introduced from 5% to 75% by 5% increments,         
                        chunkSize = NULL, # an optional integer argument giving the average length that you want          
                        # missing data gaps to be. Only used if typeMissing = "evenChunks" or "randChunks."          
                        # If a value isn't supplied, the function assumes chunk length is 5% of length of time series         
                        ...){
  if (is.null(propMiss)) {    # if "propMiss" is not provided, set it to be a vector from .05:.95 by .05    
    propMiss_f <- seq(0.05, 0.75, by = .05)  
  } else {
    propMiss_f <- propMiss  
  }    
  
  if (is.null(chunkSize)) {    
    # if "chunkSize" is not provided, set it to be 5% of length of time series (rounded to nearest integer)    
    chunkSize_f <- round(length(timeSeries)*.05)    
  } else {    
    chunkSize_f <- chunkSize  
  }    
  
  ## if you want to generate randomly missing data  
  if (typeMissing == "random") {    ## get the missing data time series for each value of missingness (returns a     
    # list where each element of the time series is increasingly missing more and more data)    
    missingDat_list <- lapply(X = propMiss_f,                               
                              FUN = function(X)
                                replace(timeSeries,                         
                                        list = sample(x = 1:length(timeSeries),
                                                      size = X*length(timeSeries),          
                                                      replace = FALSE),                                        
                                        values = NA)    )    
    
    names(missingDat_list) <- paste0("propMissAct_",  
                                     lapply(X = missingDat_list, 
                                            FUN = function(x) round(sum(is.na(x))/length(x),2)
                                     ))
  }     
  
  ## if you want to generate missing data in chunks (relatively evenly spaced)  
  if (typeMissing == "evenChunks") {
    # the 'average' chunk size is stored in "chunkSize_f"
    # iterate through the different missingness proportions (stored in "propMiss_f")
    for (i in 1:length(propMiss_f)) {
      #is the chunk size an equal or smaller amount  proportions of missing data that are input into the function
      if (chunkSize_f/length(timeSeries) > propMiss_f[i]) {
        print(paste0("there are no missing data generated for the argument resulting in ",propMiss_f[i]*100,
                     "% missing data, because a single gap results in a greater proportion of missing data than indicated for this value of the 'propMiss' argument")) 
        next      }     
      
      # calculate the shape parameter (the number of gaps that will need to
      # exist to get the specified proportion of missing data in propMiss_f[i]     
      n_events_i <- (propMiss_f[i] * length(timeSeries))/chunkSize_f      # calculate the scale parameter (what is the mean time between events that
      
      # needs to exist to get the specified proportion of missing data in propMiss_f[i]   
      beta_i <- (length(timeSeries)-(n_events_i * chunkSize_f))/n_events_i            
      
      # given the average distance between chunks calculated above, when should       
      # the next gap occur? (returns the the number of intervals needed for the number of events for this i)      
      intervals_i <- round(rgamma(n = n_events_i, shape = 1, scale = beta_i))            
      
      # make the length of the gap slightly random (draw from a poisson distribution) %%I think this is the right distribution to use?      
      chunkSize_i <- rpois(n = n_events_i, lambda = chunkSize_f)        
      
      ## now replace values in the timeSeries according in chunks of size chunkSize_i along intervals intervals_i      
      # get the indices of the starts and ends of each gap      
      indices_temp <- c(rbind(intervals_i, (chunkSize_i-1)))      
      indices_i <- c(indices_temp[1])      
      for (j in 2:length(indices_temp)) {        
        indices_i <- append(indices_i, values = indices_i[j-1] + indices_temp[j])      
      }      
      
      # if the indices are longer than the length of the timeSeries, then "recycle" them to start again in the time series      
      indices_i[indices_i > length(timeSeries)] <- indices_i[indices_i > length(timeSeries)] - length(timeSeries)      
      
      indices_seq <- NULL      
      for (k in seq(from = 2, to = length(indices_i), by = 2)) {        
        indices_seq <- append(indices_seq, indices_i[k-1]:indices_i[k])      
      }          
      
      ## store the output
      if (i == 1) {       
        missingDat_list <- list(replace(timeSeries, list = indices_seq, values = NA))     
      } else {       
        missingDat_list[[i]] <- replace(timeSeries, list = indices_seq, values = NA)     
      }          
    }    
    names(missingDat_list) <- paste0("propMissAct_",  
                                     lapply(X = missingDat_list, 
                                            FUN = function(x) round(sum(is.na(x))/length(x),2)
                                     ))
                                     }
  
  ## if you want to generate missing data in chunks (randomly spaced)  
  if (typeMissing == "randChunks") {    
     ## loop through values of "propMissing_f"
     for (i in 1:length(propMiss_f)) {
      n <- length(timeSeries)
      
      # declare transition matrix
      M <- matrix(nrow = 2, ncol = 2)
      
      # first row of transition matrix
      # determines 0 -> 0 and 0 -> 1
      M[1, ] <- c(1 - (1 / chunkSize_f), 1 / chunkSize_f)
      
      # stationary distribution
      # determined by desired proportion of missingness
      v <- c(propMiss_f[i], 1 - propMiss_f[i])
      
      # Find transition probability 1 -> 0
      # this is based on the stationary distribution
      M[2, 1] <- M[1, 2] * v[1] / v[2]
      
      # complete the matrix
      M[2, 2] <- 1 - M[2, 1]
      
      # generate series of states
      X <- vector(mode = "double", length = n)
      X[1] <- 2 # always start with an observed value
      
      for(t in 2:n){
        p_t <- M[X[t - 1], ]
        X[t] <- which(rmultinom(1, size = 1, prob = p_t) == 1)
      }
      ## is X only 0s? 
      if (sum(X != 1) == 0) {
        while (sum(X != 1) == 0) {
          # if yes, redo the process again until you get at least one '1'
          # generate series of states
          X <- vector(mode = "double", length = n)
          X[1] <- 2 # always start with an observed value
          
          for(t in 2:n){
            p_t <- M[X[t - 1], ]
            X[t] <- which(rmultinom(1, size = 1, prob = p_t) == 1)
          }
        }
      }
      
      # now change to 0s and 1s
      missingVec <- X - 1
      
      # remove the "0" values in the timeSeries string
      missingDat_temp <- replace(timeSeries, list = which(missingVec == 0), values = NA)
   
      ## store the results
      if (i == 1) {        
        missingDat_list <- list(missingDat_temp)      
      } else {        
        missingDat_list[[i]] <- missingDat_temp     
      }     
      
    }
    names(missingDat_list) <- paste0("propMissAct_",  
                                       lapply(X = missingDat_list, 
                                              FUN = function(x) round(sum(is.na(x))/length(x),2)
                                        ))
  }
  
  ## if you want to remove maximum and minimum values  
  if (typeMissing == "minMax") {    ## get the missing data time series for each 
    # value of missingness (returns a list where each element of the time series 
    # is increasingly missing more and more data)    
      
      # use a normal distribution with the same mean and sd of the timeSeries 
      # (assume it is normally distributed... ##AES change?)
      ## get the cutoff values above and below which values will be turned to NAs 
      # e.g. if 5% of data should be missing (propMiss_f == .05), then the cutoff 
      # values will be the values at .025 and 0.975 percentiles of a normal 
      # distribution with the mean and sd of the time series

      missingDat_list <- lapply(X = propMiss_f,                               
                              FUN = function(X)
                                replace(timeSeries,    
                                # get the indices of the timeSeries values above or below the cutoff values,  
                                        list = c(which(timeSeries < qnorm(X/2, mean = mean(timeSeries), sd = sd(timeSeries))),
                                                 which(timeSeries > qnorm(1-X/2, mean = mean(timeSeries), sd = sd(timeSeries)))),
                                                                               
                                        values = NA))    

      ## store the output
      if (i == 1) {       
        missingDat_list <- list(replace(timeSeries, list = indices_seq, values = NA))     
      } else {       
        missingDat_list[[i]] <- replace(timeSeries, list = indices_seq, values = NA)     
        }          

    
      names(missingDat_list) <- paste0("propMissAct_",  
                                       lapply(X = missingDat_list, 
                                              FUN = function(x) round(sum(is.na(x))/length(x),2)))
                                       }    
  return(missingDat_list)
  }

# testing -----------------------------------------------------------------
## load example data
## 1000 simulated AR(1) time series with 2 covariates to mimic the basic structure 
#of the GPP dataset. Each time series is stored in a list and named `y`. The 
# parameters used to simulate the data are stored in a sublist called `sim_params` 
# and include `phi` (the AR parameter), `beta` the linear coefficients for the covariates, 
# and `X`, the model matrix for the covariates. 
# gauss <- read_rds("./data/gauss_ar1_0miss_datasets.rds")

## 1000 simulated integer-valued time series generated using a Ricker population 
# model with Poisson-distributed demographic stochasticity. Each time series is 
# stored in a list and named `y`. The parameters used to simulate the data are 
# stored in a sublist called `sim_params` and include `r`, the intrinsic growth 
# rate, and `alpha`, the intraspecific competition coefficient.
ricker <- read_rds("./data/ricker_0miss_datasets.rds")

## testing with a subset of the ricker data
makeMissing(timeSeries = ricker[[1]]$y, typeMissing = "random", propMiss = c(.5, .05))
makeMissing(timeSeries = ricker[[1]]$y, typeMissing = "evenChunks", propMiss = c(.5, .05), chunkSize = 3)

## testing w/ all of the ricker data list
lapply(X = ricker, FUN = function(X) makeMissing(timeSeries = X$y, typeMissing = "random"))
