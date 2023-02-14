###some useful resources for time series modelling##
#https://www.analyticsvidhya.com/blog/2015/12/complete-tutorial-time-series-modeling/
#https://online.stat.psu.edu/stat501/lesson/14/14.1#


##packages##

library(dplyr)
library(ggplot2)
library(Hmisc)
library(dclone)
library(FSA)
library(FSAdata)
library(nlstools)
library(plotrix) 
library(PVAClone)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(plyr)
library(shinystan)
library(faux)
library(bayesplot)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(data.table)
library(Amelia)
library(MASS)
library(performance)



##### Simulating population dataset using RICKER model ####

Ricker<-function (Nt, r, K) 
{
  (Nt*exp(r*(1-(Nt/K))))
}


## different way of writing the Ricker model...##

Ricker2<-function (Nt, r, alpha) 
{
  (Nt*exp(r-alpha*Nt))
}


###Simulating a population dataset using Ricker model and demographic stochasticity generating from Poisson ##

simRicker_demoerror <- function(r, K, data,time) {
  
  for (t in 1:(time-1)){
    
    # Number of individuals after pop growth ###
    births <- Ricker2(data[t,], r, alpha)
    
    data[t+1,] <- rpois(1,births)  ###rpois adds demographic error around the deterministic part of function ###
  }
  
  return(data)
}


###

###Simulating a population dataset using Ricker model and demographic stochasticity generating nega. binomial ##


simRicker_nbinomerror <- function(r, K, data,time, theta) {
  
  for (t in 1:(time-1)){
    
    # Number of individuals after pop growth ###
    births <- Ricker2(data[t,], alpha, K)
    
    data[t+1,] <- rnbinom(1, size=theta, mu=births)  ###add negative binomial error ###
  }
  
  return(data)
}




### simulate data with poisson error and plot###

## r = 0.3083 , K = 308## Let's pick something similar to the Bird dataset...

r=  0.3083277   
alpha = 0.0009997187

time = 60  ##length of the times series###

results <- matrix(NA, nrow=time, ncol=1)
results[1,1] <- K  ## start w/ first time step at K###


simpop<-simRicker_demoerror(r, alpha, results, time)

simpopdf<-as.data.frame(simpop)
time<- rownames(simpopdf)
rownames(simpopdf) <- NULL
simpopts <- cbind(time, simpopdf)
simpopts$time<-as.numeric(simpopts$time)

simpop1<-ggplot(simpopts, aes(x=time, y= V1)) + 
  geom_point() + 
  geom_line()+
  xlab("")+
  ylab("Abundance")+
  ggtitle("Population simulation Ricker with Poisson error")



### simulate data with negative binomial error and plot###


## r = 0.3083 , K = 308## Let's pick something similar to the Bird dataset...

r=  0.3083   
alpha = 0.0009997187


theta= 17 ### dispersion parameter in negative binomial###

time = 60  ##length of the times series###

results <- matrix(NA, nrow=time, ncol=1)
results[1,1] <- K  ## start w/ first time step at K###

simpop2<-simRicker_nbinomerror(r, alpha, results, time, theta)

simpopdf2<-as.data.frame(simpop2)
time<- rownames(simpopdf)
rownames(simpopdf2) <- NULL
simpopts2 <- cbind(time, simpopdf2)
simpopts2$time<-as.numeric(simpopts2$time)

simpop2<-ggplot(simpopts2, aes(x=time, y= V1)) + 
  geom_point() + 
  geom_line()+
  xlab("")+
  ylab("Abundance")+
  ggtitle("Population simulation Ricker with Nbinom error")


##### play around with missing data in the simulated dataset####


## forloop that creates lists of data with increasing missingness and then solves for parameters with NLS###

simpoptry1<-as.data.frame(simpopts$V1)
colnames(simpoptry1) <- c('Broods')

missingProbs <- seq(0,1, by = .05)

# make an empty list to hold the output

simmissingData_list <- replicate(length(missingProbs),rep(NA, times = nrow(simpoptry1)) , simplify = FALSE)
simmodelOutput_list <- replicate(length(missingProbs), rep(NULL), simplify = FALSE)
for (i in 1:length(missingProbs)) {
  ## sequentially remove larger proportions of the data
  simmissingdf_i <- as.data.frame(lapply(simpoptry1,
                                      function(cc) cc[sample(c(TRUE, NA),
                                                             prob = c(1-missingProbs[i],
                                                                      missingProbs[i]),
                                                             size = length(cc),
                                                             replace = TRUE) ]))
  simmissingdf_i$popplus1 <- Lag(simmissingdf_i$Broods, -1)
  ## save the missing data in a list
  simmissingData_list[[i]] <- simmissingdf_i
  names(simmissingData_list)[i] <- paste0(missingProbs[i],"_Missing")
  ## fit an NLS model to this "i" level of missing data
  simRickermissing_i <-nls(popplus1~Broods*exp(r-alpha*Broods),start=list(r=0.1, alpha= 1/200), data=simmissingData_list[[i]])
  simmodelOutput_list[[i]] <- summary(simRickermissing_i)$coefficients
  names(simmodelOutput_list)[i] <- paste0(missingProbs[i],"_Missing")
}

####
simRickermissingdf<-plyr::ldply(simmodelOutput_list, data.frame)

paramnames<-c("r", "alpha")

paramnames1<-rep(paramnames,length(unique(simRickermissingdf$.id)))

simRickermissingdf2<-cbind(param=paramnames1,simRickermissingdf)

simRickermissingdf3<-separate(data = simRickermissingdf2, col = .id, into = c("propmissing", "right"), sep = "_")


###plot estimate and Std error of estimates over missing-ness##

ggplot(data=simRickermissingdf3, aes(x=as.numeric(propmissing), y=Estimate))+
  facet_wrap(~param, scales="free")+
  geom_errorbar(aes(ymin=Estimate-`Std..Error`, ymax=Estimate+`Std..Error`), width=.02)+
  geom_point(size=3)+
  geom_hline(data=hline_dat, aes(yintercept=value), colour="salmon")+
  theme_classic()+
  xlab("Percent of Missing Data")+
  ylab("Parameter Estimate")
