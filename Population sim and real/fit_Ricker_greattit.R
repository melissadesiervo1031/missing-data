# Run NLS and STAN model on actual Great tit population data
# Author: Melissa DeSiervo

##packages##

library(plyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(here)
library(tseries)
library(tis)

##upload the great tit population data##

titpop<- read.csv(here("data/Wytham_tits.csv"), header=T, check.names=FALSE)


##make dataset a time series object###

head(titpop)

titpopts<-ts(titpop$Broods, start=1960, end=2018, frequency=1)


###calculate lambda between years###

##add variable t+1##

titpop$poptplus1 <- Lag(titpop$Broods, k=1)

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

