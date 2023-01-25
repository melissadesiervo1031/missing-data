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

##upload the great tit population data##

titpop<- read.csv(here("data/Wytham_tits.csv"), header=T, check.names=FALSE)


##make dataset a time series object###

head(titpop)

titpopts<-ts(titpop$Broods, start=1960, end=2018, frequency=1)


##plot population time series###
tittimeseries<-ggplot(titpop, aes(x=as.numeric(as.character(Year)), y=Broods)) + 
  geom_point() + 
  geom_line()+
  xlab("")+
  ylab("Abundance")+
  ggtitle("Great Tit population")


###calculate lambda between years###

##add variable t+1##

titpop$poptplus1 <- Lag(titpop$Broods, -1)

##add variable lambda t+1/t ##

titpop2<-titpop %>% mutate(lambdayear=poptplus1/Broods)%>% mutate(percapgrowth=(poptplus1-Broods)/Broods)%>% mutate(r=log(lambdayear))

##get rid of last year##
titpop3<-complete.cases(titpop2)

titpop3<-titpop2[complete.cases(titpop2), ]


#plot population growth rate as a function of population size###

ddplotpercap<-ggplot(titpop2, aes(x=Broods, y=percapgrowth)) + 
  geom_point() + 
  geom_smooth(method=lm, se=FALSE)+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "red", size=1)+
  xlab("Population size")+
  ylab("Per capita pop growth rate")


##estimates with NLS #

Rickermodel2 <-nls(poptplus1~Broods*exp(r*(1-(Broods/K))),start=list(r=0.1, K = 200), data=titpop)

Rickermodel3 <-nls(poptplus1~Broods*exp(r-alpha*Broods),start=list(r=0.1, alpha= 1/200), data=titpop)

## r = 0.3083 , K = 308## (alpha = 0.000997)

##test whether a poisson or a neg binomial distribution is better ###

poissontit<-glm(poptplus1~Broods, family="poisson", data=titpop)

nbinomtit<-glm.nb(poptplus1~Broods, data = titpop)

check_overdispersion(poissontit)

##negative binomial is better###

###estimates with STAN###

###https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html#mechanistic-model-the-lotka-volterra-equations###
##https://www.weirdfishes.blog/blog/fitting-bayesian-models-with-stan-and-r/

# Model code for STAN -----------------------------------------------------


####
tit_dat<- list(y=titpop$Broods, N=length(titpop$Broods))

tit_fit<- stan(file="Time series/ricker3_.stan", data= tit_dat, iter = 1000, chains = 4)

print(tit_fit,digits=5)

#### r =  0.34482 , alpha =0.00112 (K = 307.9))

#how did we do?
titpred<- titpop$pop[-length(titpop$pop)]*exp(0.79*(1-(titpop$pop[-length(titpop$pop)]/211)))
diff_titpred<- titpred-titpop$pop[-length(titpop$pop)]
diff_titpop <- titpop$pop[-1]-titpop$pop[-length(titpop$pop)]

plot(diff_titpop,diff_titpred)

#### Matches the NLS estimates###

#### calculate estimates with maximum likilihood approach (fitting_ricker_model.R)

#years##

y=titpop$Broods

years=max(titpop$Year)-min(titpop$Year)+1

# provide X matrix
X <- cbind(
  rep(1, years),
  y
)


fitLL <- optim(par = c(1, 0), fn = ricker_Pois_neg_ll, y = y, X = X, hessian = T)

##not exactly matching with NLS OR STAN, ## r =  0.491790068 , alpha = 0.001561374 (K = 315)




#### solve for estimates with ARIMA ###

##helpful links https://ionides.github.io/531w16/final_project/Project19/final.html ##
#https://towardsdatascience.com/state-space-model-and-kalman-filter-for-time-series-prediction-basic-structural-dynamic-linear-2421d7b49fa6#
#https://lbelzile.github.io/timeseRies/state-space-models-and-the-kalman-filter.html



########play around with missing data#############

titpopdf1<-as.data.frame(titpop$Broods)
colnames(titpopdf1) <- c('Broods')


###take away 15% of estimates##

missingdf1<-as.data.frame(lapply(titpopdf1, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.85, 0.15), size = length(cc), replace = TRUE) ]))

missingdf11<-missingdf1 %>% mutate(popplus1=Lag(titpop$Broods, -1))


#####
## NLS can handle missing data, but will accumulate bias##

Rickermodel22 <-nls(popplus1~Broods*exp(r*(1-(Broods/(1/a)))),start=list(r=0.1, a= 1/200), data=missingdf11)



