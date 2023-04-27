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
      
      # make the length of the gap slightly random (draw from a poisson distribution) %%I think this is the right distribution to use?      
      # get a vector of possible gap sizes, then pick the first 
      chunkSize_TEMP <- rpois(n = 100, lambda = chunkSize_f)        
      # make sure that chunks are at least 1  unit long
      chunkSize_TEMP[chunkSize_TEMP==0] <- 1
      # get the correct number of chunks that correspond to the amount of missing 
      # data that the propMiss_f[i] argument specificies
      gapsNeeded <- propMiss_f[i]*length(timeSeries)
      
      chunkMat <- matrix(NA, nrow = 100, ncol = 100)
      for (k in 1:length(chunkSize_TEMP)) {
        chunkMat[,k] <- c(chunkSize_TEMP[1:k], rep(NA, length.out = 100-k))
      }
      # find out how many of the chunks we need to get the amount of missing data we want
      chunkSize_i <- chunkSize_TEMP[1:which.min(abs(colSums(chunkMat, na.rm = TRUE) - gapsNeeded))]
    
      # calculate the shape parameter (the number of gaps that will need to
      # exist to get the specified proportion of missing data in propMiss_f[i]     
      n_events_i <- length(chunkSize_i)  
      
      # calculate the scale parameter (what is the mean time between events that
      # needs to exist to get the specified proportion of missing data in propMiss_f[i]   
      beta_i <- round((length(timeSeries) - sum(chunkSize_i))/n_events_i)         
      
      # given the average distance between chunks calculated above, when should       
      # the next gap occur? (returns the the number of intervals needed for the number of events for this i)      
      # intervals are random, but make it so that the length of the events/iterations can't be longer than the length of the time series
      repeat {
        intervals_i <- round(rgamma(n = n_events_i, shape = 1, scale = beta_i))  
        if (sum(intervals_i, chunkSize_i) <= length(timeSeries)) {
          break
        }
      }          
      
      ## now replace values in the timeSeries according in chunks of size chunkSize_i along intervals intervals_i      
      # get the indices of the starts and ends of each gap      
      indices_temp <- c(rbind((intervals_i+1), (chunkSize_i-1)))      
      indices_j <- c(indices_temp[1])      
      for (j in 2:length(indices_temp)) {        
        indices_j <- append(indices_j, values = indices_j[j-1] + indices_temp[j])      
      }      
    
      indices_seq <- NULL      
      for (k in seq(from = 2, to = length(indices_j), by = 2)) {        
        indices_seq <- append(indices_seq, indices_j[k-1]:indices_j[k])      
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
ricker <- readRDS("./data/ricker_0miss_datasets.rds")

## testing with a subset of the ricker data
makeMissing(timeSeries = ricker[[1]]$y, typeMissing = "random", propMiss = c(.5, .05))
makeMissing(timeSeries = ricker[[1]]$y, typeMissing = "evenChunks", propMiss = c(.5, .05), chunkSize = 3)

## testing w/ all of the ricker data list
lapply(X = ricker, FUN = function(X) makeMissing(timeSeries = X$y, typeMissing = "random"))


# Testing for MCAR vs MNAR ------------------------------------------------
library(ggpubr)

# use the simulated gaussian AR1 data
gausSim <- readRDS("./data/gauss_ar1_0miss_datasets.rds")
# for now, get just the first dataset
gausSim <- data.frame("y" = gausSim[[1]]$y, 
                      "time" = 1:length(gausSim[[1]]$y))

# get 40% missing completely at random 
gausSim$y_randMiss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "random", propMiss = .4)))

# get 40% missing in evenly spaced chunks
gausSim$y_evenChunkMiss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "evenChunks", propMiss = .4, chunkSize = 5)))

# get 40% missing in randomly spaced chunks
gausSim$y_randChunkMiss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "randChunks", propMiss = .4, chunkSize = 5)))

# get 40% missing in minmax of data
gausSim$y_minMaxMiss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "minMax", propMiss = .4)))

# make histograms like Dusty made (looking at histogram of values in the missing
# data time series vs. expected distribution from the AR1 process)
# no missing data
noMiss_line <- ggplot(data = gausSim, aes(x = time, y = y)) + 
  geom_line() + 
  geom_point(size = 1) +
  ggtitle("No Missing Data") + 
  theme_classic()

noMiss_hist <- ggplot() + 
  geom_histogram(data = gausSim, aes(y, after_stat(density)), fill = "grey", color = "darkgrey") + 
  geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)), 
                x = seq(-5,7,.1)), color = "blue") + 
  theme_classic()

# randomly missing
randMiss_line <- ggplot(data = gausSim, aes(x = time, y = y_randMiss)) + 
  geom_line() + 
  geom_point(size = 1) +
  ggtitle("Missing Completely At Random") + 
  theme_classic()

randMiss_hist <- ggplot() + 
  geom_histogram(data = gausSim, aes(y_randMiss, after_stat(density)), fill = "grey", color = "darkgrey") + 
  geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)), 
                x = seq(-5,7,.1)), color = "blue") + 
  theme_classic()

# missing evenly spaced chunks
evenChunkMiss_line <- ggplot(data = gausSim, aes(x = time, y = y_evenChunkMiss)) + 
  geom_line() + 
  geom_point(size = 1) +
  ggtitle("Missing Even Chunks") + 
  theme_classic()

evenChunkMiss_hist <- ggplot() + 
  geom_histogram(data = gausSim, aes(y_evenChunkMiss, after_stat(density)), fill = "grey", color = "darkgrey") + 
  geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)), 
                x = seq(-5,7,.1)), color = "blue") + 
  theme_classic()

# missing randomly spaced chunks
randChunkMiss_line <- ggplot(data = gausSim, aes(x = time, y = y_randChunkMiss)) + 
  geom_line() + 
  geom_point(size = 1) +
  ggtitle("Missing Randomly Spaced Chunks") + 
  theme_classic()

randChunkMiss_hist <- ggplot() + 
  geom_histogram(data = gausSim, aes(y_randChunkMiss, after_stat(density)), fill = "grey", color = "darkgrey") + 
  geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)), 
                x = seq(-5,7,.1)), color = "blue") + 
  theme_classic()

# missing at min and max
minMaxMiss_line <- ggplot(data = gausSim, aes(x = time, y = y_minMaxMiss)) + 
  geom_line() + 
  geom_point(size = 1) +
  ggtitle("Missing Min and Max") + 
  ylim(c(-2.5,5)) +
  theme_classic()

minMaxMiss_hist <- ggplot() + 
  geom_histogram(data = gausSim, aes(y_minMaxMiss, after_stat(density)), fill = "grey", color = "darkgrey") + 
  geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)), 
                x = seq(-5,7,.1)), color = "blue") + 
  theme_classic()

ggarrange(noMiss_line, noMiss_hist, 
          randMiss_line, randMiss_hist, 
          evenChunkMiss_line, evenChunkMiss_hist,
          randChunkMiss_line, randChunkMiss_hist,
          minMaxMiss_line, minMaxMiss_hist,
          widths = c(.75,.25, .75,.25, .75,.25, .75,.25, .75,.25), 
          ncol = 2, 
          nrow = 5)
