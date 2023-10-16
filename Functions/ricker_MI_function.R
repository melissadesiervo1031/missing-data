########################################################################
# Function will fill in missing data with multiple imputation using Amelia
# and then fit the Ricker Model


# Load packages
library(here)
library(stats)
#library(forecast)
#library(tidyverse)
#library(lubridate)
library(Amelia)




#methods: forward, backward, dual, averaging
#' Title Fit Ricker Model with Multiple Imputation using Amelia
#'
#' @param y vector of population time series data
#' @param imputationsnum number of imputations to average together, more will take longer but have greater accuracy 
#' @param fam either fit the models using poisson ("poisson") or negative binomial ("neg_binom")
#' @param method Each missing data point will be filled in twice by Amelia, once as y_t and once as y_t-1, should Amelia only predict y_t ("forward"), only predict y_t-1 ("backward"), use both and allow discrepancies ("dual"), or use both and average out discrepancies ("averaging")
#'
#' @return
#' @export
#'
#' @examples
fit_ricker_MI<-function(y, imputationsnum, fam = "poisson", method="dual"){
  
  # get length of time series
  n <- length(y)
  
  # make data frame in prep for multiple imputation
  simmissingdf=cbind.data.frame(1:(n-1),y[2:n],y[1:(n - 1)])
  names(simmissingdf)<-c("time","yt","yt1")
  
  # do the multiple imputation 
  # include lags? (doesn't fix missing data)
  # how to include longer stretches of NA's
  # is it a problem that we are basically using Amelia to predict the same value twice?
  
  # get rid of NAs both forwards and backwards, without averaging
  if(method=="dual"){
    amelia1sim<-amelia(simmissingdf, m=imputationsnum, ts="time")
    for(i in 1:imputationsnum){
      while(length(which(is.na(amelia1sim$imputations[[i]]$yt)))>0){
        # get where the NAs are
        navec<-which(is.na(amelia1sim$imputations[[i]]$yt))
        # replace some NAs with values that should be the same
        amelia1sim$imputations[[i]]$yt[navec]<-amelia1sim$imputations[[i]]$yt1[navec+1]
        amelia1sim$imputations[[i]]$yt1[navec]<-amelia1sim$imputations[[i]]$yt[navec-1]
        # do another amelia round to fill in any more NAs if possible
        am2<-amelia(amelia1sim$imputations[[i]],m=1,ts="time")
        amelia1sim$imputations[[i]]<-am2$imputations$imp1
      }
    }
  }
  
  if(method=="forward"){
    amelia1sim<-amelia(simmissingdf, m=imputationsnum, ts="time")
    for(i in 1:imputationsnum){
      #correct it to only predict forward, that is only predict yt, and then fill the yt1
      fill1=which(amelia1sim$imputations[[i]][1:(n-2),2]!=amelia1sim$imputations[[i]][2:(n-1),3])
      amelia1sim$imputations[[i]][fill1+1,3]=amelia1sim$imputations[[i]][fill1,2]
      
      while(length(which(is.na(amelia1sim$imputations[[i]]$yt)))>0){
        # get where the NAs are
        navec<-which(is.na(amelia1sim$imputations[[i]]$yt))
        # replace some NAs with values that should be the same
        #amelia1sim$imputations[[i]]$yt[navec]<-amelia1sim$imputations[[i]]$yt1[navec+1]
        amelia1sim$imputations[[i]]$yt1[navec]<-amelia1sim$imputations[[i]]$yt[navec-1]
        # do another amelia round to fill in any more NAs if possible
        am2<-amelia(amelia1sim$imputations[[i]],m=1,ts="time")
        amelia1sim$imputations[[i]]<-am2$imputations$imp1
        #correct it to only predict forward, that is only predict yt, and then fill the yt1
        fill1=which(amelia1sim$imputations[[i]][1:(n-2),2]!=amelia1sim$imputations[[i]][2:(n-1),3])
        amelia1sim$imputations[[i]][fill1+1,3]=amelia1sim$imputations[[i]][fill1,2]
        
      }
    }
  }
  
  if(method=="backward"){
    amelia1sim<-amelia(simmissingdf, m=imputationsnum, ts="time")
    for(i in 1:imputationsnum){
      #correct it to only predict backward, that is only predict yt1, and then fill the yt
      fill1=which(amelia1sim$imputations[[i]][1:(n-2),2]!=amelia1sim$imputations[[i]][2:(n-1),3])
      amelia1sim$imputations[[i]][fill1,2]=amelia1sim$imputations[[i]][fill1+1,3]
      
      while(length(which(is.na(amelia1sim$imputations[[i]]$yt)))>0){
        # get where the NAs are
        navec<-which(is.na(amelia1sim$imputations[[i]]$yt))
        # replace some NAs with values that should be the same
        amelia1sim$imputations[[i]]$yt[navec]<-amelia1sim$imputations[[i]]$yt1[navec+1]
        #amelia1sim$imputations[[i]]$yt1[navec]<-amelia1sim$imputations[[i]]$yt[navec-1]
        # do another amelia round to fill in any more NAs if possible
        am2<-amelia(amelia1sim$imputations[[i]],m=1,ts="time")
        amelia1sim$imputations[[i]]<-am2$imputations$imp1
        #correct it to only predict backward, that is only predict yt1, and then fill the yt
        fill1=which(amelia1sim$imputations[[i]][1:(n-2),2]!=amelia1sim$imputations[[i]][2:(n-1),3])
        amelia1sim$imputations[[i]][fill1,2]=amelia1sim$imputations[[i]][fill1+1,3]
        
      }
    }
  }
  
  if(method=="averaging"){
    amelia1sim<-amelia(simmissingdf, m=imputationsnum, ts="time")
    for(i in 1:imputationsnum){
     
      while(length(which(is.na(amelia1sim$imputations[[i]]$yt)))>0){
        # get where the NAs are
        navec<-which(is.na(amelia1sim$imputations[[i]]$yt))
        # replace some NAs with values that should be the same
        amelia1sim$imputations[[i]]$yt[navec]<-amelia1sim$imputations[[i]]$yt1[navec+1]
        amelia1sim$imputations[[i]]$yt1[navec]<-amelia1sim$imputations[[i]]$yt[navec-1]
        # do another amelia round to fill in any more NAs if possible
        am2<-amelia(amelia1sim$imputations[[i]],m=1,ts="time")
        amelia1sim$imputations[[i]]<-am2$imputations$imp1
      }
      #correct it to average the predicted yt and yt1 that don't match
      fill1=which(amelia1sim$imputations[[i]][1:(n-2),2]!=amelia1sim$imputations[[i]][2:(n-1),3])
      fill2=(amelia1sim$imputations[[i]][fill1,2]+amelia1sim$imputations[[i]][fill1+1,3])/2
      amelia1sim$imputations[[i]][fill1,2]=fill2
      amelia1sim$imputations[[i]][fill1+1,3]=fill2
    }
  }
  
  
  
  fit=list()
  # fit model over all imputations
  for(i in 1:imputationsnum){
    
    # compile into sliced dataframe in preparation for 
    dat <- data.frame(
      yt = amelia1sim$imputations[[i]][2:n,2],
      ytm1 = amelia1sim$imputations[[i]][1:(n - 1),2]
    )
    
    # fit ricker model with poisson
    if(fam == "poisson"){
      fit[[i]] <- glm(
        round(yt) ~ round(ytm1), data = dat, 
        family = poisson,
        offset = log(ytm1)
      )
    }
    
    # or fit ricker model with negbinom
    if(fam == "neg_binom"){
      fit[[i]] <- MASS::glm.nb(
        yt ~ ytm1 + offset(log(ytm1)), 
        data = dat
      )
    }
    
  }
  
  
  # Averages models into 1 result and simplifies
  estims1=sapply(fit,coef,simplify=T)
  estims=rowMeans(estims1)
  
  cis1=sapply(fit,confint,simplify=T)
  cis=rowMeans(cis1)
  
  ses1=sapply(fit,function(X) sqrt(diag(vcov(X))),simplify=T)
  ses=rowMeans(ses1)
  
  # compile objects and rename
  #estims <- coef(fit)
  names(estims) <- c("r", "alpha")
  #cis <- confint(fit)
  names(cis) <- rep(c("r", "alpha"),2)
  #ses <- sqrt(diag(vcov(fit)))
  names(ses) <- names(estims)
  
  # return as a list
  return(list(
    estim = estims * c(1, -1),
    se = ses,
    lower = as.double(cis[c(1,4)] * c(1, -1)),
    upper = as.double(
      cis[c(3,2)] * c(1, -1)
    )
  ))
  
}

#example
poiss_sims <- readRDS("/Users/amypatterson/Documents/Laramie_postdoc/Missing data TS/missing-data/data/missingDatasets/pois_sim_randMiss_A.rds")
y=poiss_sims[[1]]$y[[1]]
imputationsnum=5

ans1=fit_ricker_MI(y,imputationsnum = 5,method="dual")
ans2=fit_ricker_MI(y,imputationsnum = 5,method="forward")
ans3=fit_ricker_MI(y,imputationsnum = 5,method="backward")
ans4=fit_ricker_MI(y,imputationsnum = 5,method="averaging")


