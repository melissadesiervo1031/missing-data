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
                        # generated ("random" (with varying degrees of autocorrelation) or "minMax")         
                        propMiss = NULL, # an optional vector argument that gives the proportion(s)          
                        # of the data that should be missing; if not specified, then increasingly          
                        # more missingness is introduced from 5% to 75% by 5% increments,         
                        autoCorr= NULL, # an optional argument between 0 and 1 (non-inclusive) giving the degree of 
                        # autocorrelation in missingness that the timeseries will have. (0 is no autocorrelation, .99 is maximum autocorrelation)
                        # Only used if typeMissing = "random". If a value isn't supplied, the function assumes no autocorrelation (0)        
                        ...){
  if (is.null(propMiss)) {    # if "propMiss" is not provided, set it to be a vector from .05:.75 by .05    
    propMiss_f <- seq(0.05, 0.75, by = .05)  
  } else {
    propMiss_f <- propMiss  
  }    
  
  if (is.null(autoCorr)) {    
    # if "chunkSize" is not provided, set it to be 5% of length of time series (rounded to nearest integer)    
    autoCorr_f <- .5
  } else {    
    autoCorr_f <- autoCorr  
  }    
  
  ## if you want to generate missing data in chunks (randomly spaced)  
  if (typeMissing == "random") {      
    ## loop through values of "propMissing_f"
    for (i in 1:length(propMiss_f)) {
      
      n <- length(timeSeries)
      
      # declare transition matrix
      M <- matrix(nrow = 2, ncol = 2)
      
      # populate transition matrix based on Gharib et al. 2014 ("Characterization of Markov-Bernoulli geometric distribution related to random sums")
      #         X_i+1
      #          0                 1
      #          ______________________
      #       0| 1-(1-rho)p      (1-rho)p
      # X_i   1| (1-rho)(1-p)    rho + (1-rho)p
      #
      # with the initial distribution P(X_1=1) = p = 1-P(X_1 = 0)
      # (so p = 1-propMiss_f[i])
      p <- 1-propMiss_f[i]
      
      ## rho is the correlation coefficient between values at times t and t+1
      # 0 indicates maximum negative autocorrelation
      # 1 indicates the maximum positive autocorrelation
      # 0.5 indicates no autocorrelation
      
      # first row of transition matrix
      # determines 0 -> 0 and 0 -> 1
      M[1, ] <- c(1 - (1 - autoCorr_f)*p, (1 - autoCorr_f)*p)
      
      # find second row of matrix
      # determines 1 -> 0 and 1 -> 1
      M[2,] <- c((1 - autoCorr_f)*(1-p), autoCorr_f + (1 - autoCorr_f)*p)
      
      # check that there are no negative or >1 probabilities in the transition matrix 
      # (i.e. not a feasible combination of chunk size and desired proportion of missingness)
      if (sum(M > 1) > 0 | sum(M < 0) > 0) {
        missingDat_temp <- c("no viable transition matrix for this combination of autocorrelation and proportion missing")
        
        ## store the results
        if (i == 1) {        
          missingDat_list <- list(missingDat_temp)      
        } else {        
          missingDat_list[[i]] <- missingDat_temp     
        }     
        next
      }
      
      # generate series of states
      X <- vector(mode = "double", length = n)
      X[1] <- 2 # always start with an observed value
      
      for(t in 2:n){
        p_t <- M[X[t - 1], ]
        X[t] <- which(rmultinom(1, size = 1, prob = p_t) == 1)
      }
      ## is X only 1s? 
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
      ## is X only 2s? 
      if (sum(X != 2) == 0) {
        while (sum(X != 2) == 0) {
          # if yes, redo the process again until you get at least one '2'
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
      
      # ## calculate the actual autocorrelation of the time series in this iteration (with a time lag of 1)
      actualAutoCorr <- acf(x = missingVec, plot = FALSE, lag.max = 1)$acf[,,1][2]
      
      actualAutoCorr_while <- actualAutoCorr
      while (actualAutoCorr_while < 0) {
        # generate series of states
        X <- vector(mode = "double", length = n)
        X[1] <- 2 # always start with an observed value
        
        for(t in 2:n){
          p_t <- M[X[t - 1], ]
          X[t] <- which(rmultinom(1, size = 1, prob = p_t) == 1)
        }
        ## is X only 1s? 
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
        ## is X only 2s? 
        if (sum(X != 2) == 0) {
          while (sum(X != 2) == 0) {
            # if yes, redo the process again until you get at least one '2'
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
        
        # ## calculate the actual autocorrelation of the time series in this iteration (with a time lag of 1)
        actualAutoCorr_while <- acf(x = missingVec, plot = FALSE, lag.max = 1)$acf[,,1][2]
        
      }
      
      # ## quick checks:
      # # average size of 'chunks' of missing data
      # ones <- which(missingVec==1)
      # mean(diff(ones)[diff(ones) != 1]-1)
      # # proportion of missing data
      # table(missingVec)[1]/n
      # # stationary distribution
      # M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M %*% M
      # 
      # remove the "0" values in the timeSeries string
      missingDat_temp <- replace(timeSeries, list = which(missingVec == 0), values = NA)
      
      ## store the results
      if (i == 1) {        
        missingDat_list <- list(missingDat_temp) 
      } else {        
        missingDat_list[[i]] <- missingDat_temp     
      }     
      # add a name w/ the true proportion of missing and the true autocorrelation value for a lag of 1
      names(missingDat_list)[i] <- paste0("propMissAct_",  
                                          round(sum(is.na(missingDat_list[[i]]))/length(missingDat_list[[i]]),2),
                                          "_autoCorr_",
                                          round(acf(x = missingVec, plot = FALSE, lag.max = 1)$acf[,,1][2],2)
      )
      
    }
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
    
    
    names(missingDat_list) <- paste0("propMissAct_",  
                                     lapply(X = missingDat_list, 
                                            FUN = function(x) round(sum(is.na(x))/length(x),2)))
  }    
  
  return(missingDat_list)
}

# testing for functionality -----------------------------------------------------------------
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
# ricker <- readRDS("./data/ricker_0miss_datasets.rds")
# 
# ## testing with a subset of the ricker data
# # for missing at random (chunk size = 1)
# makeMissing(timeSeries = ricker[[1]]$y, typeMissing = "random", propMiss = c(.5, .05), autoCorr = .003)
# # for missing at random but autocorrelated (chunk size = 3)
# makeMissing(timeSeries = ricker[[1]]$y, typeMissing = "random", propMiss = c(.5, .05))
# 
# ## testing w/ all of the ricker data list
# lapply(X = ricker, FUN = function(X) makeMissing(timeSeries = X$y, typeMissing = "random", chunkSize = 1))


# # Testing for MCAR vs MNAR ------------------------------------------------
# library(ggpubr)
# 
# # use the simulated gaussian AR1 data
# gausSim <- readRDS("./data/gauss_ar1_0miss_datasets.rds")
# # for now, get just the first dataset
# gausSim <- data.frame("y" = gausSim[[1]]$y, 
#                       "time" = 1:length(gausSim[[1]]$y))
# 
# # get 40% missing completely at random 
# gausSim$y_randMiss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "random", propMiss = .4, autoCorr = .01)))
# 
# # get 40% missing in highly autocorrelated chunks
# gausSim$y_randChunkMiss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "random", propMiss = .4, autoCorr = .9)))
# 
# # get 40% missing in minmax of data
# gausSim$y_minMaxMiss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "minMax", propMiss = .4)))
# 
# # make histograms like Dusty made (looking at histogram of values in the missing
# # data time series vs. expected distribution from the AR1 process)
# # no missing data
# noMiss_line <- ggplot(data = gausSim, aes(x = time, y = y)) + 
#   geom_line() + 
#   geom_point(size = 1) +
#   ggtitle("No Missing Data") + 
#   theme_classic()
# 
# noMiss_hist <- ggplot() + 
#   geom_histogram(data = gausSim, aes(y, after_stat(density)), fill = "grey", color = "darkgrey") + 
#   geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)), 
#                 x = seq(-5,7,.1)), color = "blue") + 
#   theme_classic()
# 
# # randomly missing
# randMiss_line <- ggplot(data = gausSim, aes(x = time, y = y_randMiss)) + 
#   geom_line() + 
#   geom_point(size = 1) +
#   ggtitle("Missing at Random: 40% missing, autoCorr = .01") + 
#   theme_classic()
# 
# randMiss_hist <- ggplot() + 
#   geom_histogram(data = gausSim, aes(y_randMiss, after_stat(density)), fill = "grey", color = "darkgrey") + 
#   geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)), 
#                 x = seq(-5,7,.1)), color = "blue") + 
#   theme_classic()
# 
# # missing randomly spaced chunks
# randChunkMiss_line <- ggplot(data = gausSim, aes(x = time, y = y_randChunkMiss)) + 
#   geom_line() + 
#   geom_point(size = 1) +
#   ggtitle("Missing at Random: 40% missing, autoCorr = .9") + 
#   theme_classic()
# 
# randChunkMiss_hist <- ggplot() + 
#   geom_histogram(data = gausSim, aes(y_randChunkMiss, after_stat(density)), fill = "grey", color = "darkgrey") + 
#   geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)), 
#                 x = seq(-5,7,.1)), color = "blue") + 
#   theme_classic()
# 
# # missing at min and max
# minMaxMiss_line <- ggplot(data = gausSim, aes(x = time, y = y_minMaxMiss)) + 
#   geom_line() + 
#   geom_point(size = 1) +
#   ggtitle("Missing Min and Max: 40% missing") + 
#   ylim(c(-2.5,5)) +
#   theme_classic()
# 
# minMaxMiss_hist <- ggplot() + 
#   geom_histogram(data = gausSim, aes(y_minMaxMiss, after_stat(density)), fill = "grey", color = "darkgrey") + 
#   geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)), 
#                 x = seq(-5,7,.1)), color = "blue") + 
#   theme_classic()
# 
# 
# ## save the figure to file
# png(filename = "./figures/CompareMissingnessTypes_fig.png", width = 700, height = 600)
# ggarrange(noMiss_line, noMiss_hist, 
#           randMiss_line, randMiss_hist, 
#           randChunkMiss_line, randChunkMiss_hist,
#           minMaxMiss_line, minMaxMiss_hist,
#           widths = c(.75,.25, .75,.25, .75,.25, .75,.25, .75,.25), 
#           ncol = 2, 
#           nrow = 4)
# dev.off()
