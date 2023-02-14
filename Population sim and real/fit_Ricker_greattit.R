# Run NLS and STAN model on actual Great tit population data
# Author: Melissa DeSiervo

##packages##

library(dplyr)
library(ggplot2)
library(tidyverse)
library(here)
library(tseries)
library(tis)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(plyr)
library(shinystan)

##upload the great tit population data##

titpop <- read_csv('data/Wytham_tits.csv')

titpop<-as.data.frame(titpop)

##make dataset a time series object###

head(titpop)

titpopts<-ts(titpop$Broods, start=1960, end=2018, frequency=1)


###calculate lambda between years###

##add variable t+1##

titpop$poptplus1 <- c(titpop$Broods[-1], NA)

##add variable lambda t+1/t ##

titpop2<-titpop %>% mutate(lambdayear=poptplus1/Broods)%>% mutate(percapgrowth=(poptplus1-Broods)/Broods)%>% mutate(r=log(lambdayear))

##get rid of last year##
titpop3<-complete.cases(titpop2)

titpop3<-titpop2[complete.cases(titpop2), ]


##plot population time series###
tittimeseries<-ggplot(titpop, aes(x=as.numeric(as.character(Year)), y=Broods)) + 
  geom_point() + 
  geom_line()+
  xlab("")+
  ylab("Abundance")+
  ggtitle("Great Tit population")


#plot population growth rate as a function of population size###

ddplotpercap<-ggplot(titpop2, aes(x=Broods, y=percapgrowth)) + 
  geom_point() + 
  geom_smooth(method=lm, se=FALSE)+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "red", size=1)+
  xlab("Population size")+
  ylab("Per capita pop growth rate")


##Ricker model estimates with NLS # note that NLS assumes normally distributed data##
#this dataset is actually best fit with a negative binomial## 

Rickermodel3 <-nls(poptplus1~Broods*exp(r-alpha*Broods),start=list(r=0.1, alpha= 1/200), data=titpop)

summary(Rickermodel3)

##r = 0.308, alpha = 0.0009 ##

### Run STAN model for estimating coefficients# Ricker model w/ Poisson error ###

tit_dat<- list(y=titpop$Broods, N=length(titpop$Broods))

tit_fit<- stan(file="Population sim and real/ricker3_.stan", data= tit_dat, iter = 4000, chains = 4)

print(tit_fit,digits=5)

##r = 0.3474, alpha = 0.00112 ##

# examine model outputs
traceplot(tit_fit, pars=c("rmax", "alpha"))
pairs(tit_fit, pars=c("rmax", "alpha","lp__"))
stan_dens(tit_fit, pars=c("rmax", "alpha"))

fit_extract <- rstan::extract(tit_fit)


# look at posterior predictions:
get_pp <- function(fit){
  fit_extract <- rstan::extract(tit_fit)
  pp <- t(apply(fit_extract$pop_rep, 2, 
                function(x) quantile(x, probs = c(0.025, 0.5, 0.975), names = F))) %>%
    data.frame() 
  names(pp) <- c('pop_rep.lower', 'pop_rep', 'pop_rep.upper')
  
  return(pp)
  
}
