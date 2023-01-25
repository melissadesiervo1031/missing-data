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





###Simulating a population dataset using Ricker model and demographic stochasticity generating from Poisson ##

simRicker_demoerror <- function(r, K, data,time) {
  
  for (t in 1:(time-1)){
    
    # Number of individuals after pop growth ###
    births <- Ricker(data[t,], r, K)
    
    data[t+1,] <- rpois(1,births)  ###rpois adds demographic error around the deterministic part of function ###
  }
  
  return(data)
}


###

###Simulating a population dataset using Ricker model and demographic stochasticity generating nega. binomial ##


simRicker_nbinomerror <- function(r, K, data,time, theta) {
  
  for (t in 1:(time-1)){
    
    # Number of individuals after pop growth ###
    births <- Ricker(data[t,], r, K)
    
    data[t+1,] <- rnbinom(1, size=theta, mu=births)  ###add negative binomial error ###
  }
  
  return(data)
}




### simulate data with poisson error and plot###

## r = 0.3083 , K = 308## Let's pick something similar to the Bird dataset...

r=  0.3083   
K = 308

time = 60  ##length of the times series###

results <- matrix(NA, nrow=time, ncol=1)
results[1,1] <- K  ## start w/ first time step at K###


simpop<-simRicker_demoerror(r, K, results, time)

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
K = 308


theta= 17 ### dispersion parameter in negative binomial###

time = 60  ##length of the times series###

results <- matrix(NA, nrow=time, ncol=1)
results[1,1] <- K  ## start w/ first time step at K###

simpop2<-simRicker_nbinomerror(r, K, results, time, theta)

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



### estimate paramaters with our simulated dataset##

####real values   ## r = 0.3083 , K = 308##

##add variable t+1##

simpopts$poptplus1 <- Lag(simpopts$V1, -1)

fit <-nls(poptplus1~V1*exp(r*(1-(V1/K))),start=list(r=0.1, K = 200), data=simpopts)



