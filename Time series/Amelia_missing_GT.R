################### With the Wytham woods dataset################
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
library(here)
library(tseries)

#helpful resource for MI##
#https://thomasleeper.com/Rcourse/Tutorials/mi.html#


##upload the great tit population data##

titpop<- read.csv(here("data/Wytham_tits.csv"), header=T, check.names=FALSE)


##make dataset a time series object###

head(titpop)

titpopts<-ts(titpop$Broods, start=1960, end=2018, frequency=1)

titpopdf1<-as.data.frame(titpop$Broods)

colnames(titpopdf1) <- c('Broods')

###make a datsaset with 20 % missing data to play around with AMELIA##

missingdf20<-as.data.frame(lapply(titpopdf1, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.80, 0.20), size = length(cc), replace = TRUE) ]))

missingdf20_2 <- tibble::rownames_to_column(missingdf20, "Time")

missingdf20_2$Time<-as.numeric(missingdf20_2$Time)

missingdf20_3<-missingdf20 %>% mutate(popplus1=Lag(missingdf20$Broods, -1))

###practice w/ AMELIA##

amelia1<-amelia(missingdf20_2, m=5, ts="Time", polytime=1) ##5 imputed datasets###  
##with polytime = 1, considers linear effects of time# 
##might also want to play around with lags and leads###
## also play around with priors ###

plot(amelia1) ##plots of observed and imputed values##

overimpute(amelia1, var=2, main = "Observed versus Imputed Values") ##check on imputation, graphs estimates of observations vs. true values##

missmap(amelia1) ##maps out where the missing data are in time series###


##to pull out one imputation of amelia ##

amelia1[["imputations"]][["imp1"]][["Time"]]

####

## for loop to do an NLS model from Great Tit for AMELIA imputation###

# make empty lists to hold the output

amelia_list <- replicate(length(5),rep(NA, times = nrow(titpopdf1)) , simplify = FALSE)
amelia_list_t1 <- replicate(length(5),rep(NA, times = nrow(titpopdf1)) , simplify = FALSE)

df1<-data.frame(matrix(NA, nrow=nrow(titpopdf1), ncol=2))

listofdf<-list(df1, df1, df1, df1, df1)

modelOutput_list <- replicate(length(5), rep(NULL), simplify = FALSE)

###

for(i in 1:5){
  amelia_list [[i]] <-amelia1[["imputations"]][[i]]$Broods
  amelia_list_t1[[i]]<- Lag(amelia_list [[i]],-1)
  listofdf[[i]]<-data.frame(Broods=unlist(amelia_list [[i]]),popplus1=unlist(amelia_list_t1[[i]]))
  
  ## fit an NLS model to each AMELIA dataset ##
  Rickermissing_i <-nls(popplus1~Broods*exp(r-alpha*Broods),start=list(r=0.1, alpha= 1/200), data=listofdf[[i]])
  modelOutput_list[[i]] <- summary(Rickermissing_i)$coefficients
 }


#### another way to do this part using an apply function###

#lm.amelia.out <- lapply(imp.amelia$imputations, function(i) lm(y ~ x1 + x2, 
                                                               #data = i))
##model output list into a dataframe to average..##


coefdf<-do.call(rbind, modelOutput_list)
params <- rownames(coefdf)
coefdf1<-cbind(params, as.data.frame(coefdf[,1:2]))
rownames(coefdf1)<-NULL

est1<-coefdf1[,1:2]
se1<-coefdf1[,c(1,3)]

est2<-pivot_wider(est1)
##need to do model averaging across AMELIA impuations##

## Uses Rubin's rules for combining a set of results from multiply imputed datasets to reflect the average result, #
mi.meld(coefs.amelia, ses.amelia)

####


