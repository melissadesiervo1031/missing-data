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

##upload the data##

titpop <- read.csv("C:/data/Wytham_tits.csv", header=T, check.names=FALSE)



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

## r = 0.3083 , K = 308##

##test whether a poisson or a neg binomial distribution is better with data###


poissontit<-glm(poptplus1~Broods, family="poisson", data=titpop)

nbinomtit<-glm.nb(poptplus1~Broods, data = titpop)

check_overdispersion(poissontit)

##negative binomial is better###



###estimates with STAN###

###https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html#mechanistic-model-the-lotka-volterra-equations###
##https://www.weirdfishes.blog/blog/fitting-bayesian-models-with-stan-and-r/

# Model code for STAN -----------------------------------------------------


####
setwd("C:Population_great tit")
tit_dat<- list(y=titpop$Broods, N=length(titpop$Broods))

tit_fit<- stan(file="ricker3_.stan", data= tit_dat, iter = 1000, chains = 4)

print(tit_fit)

###r = 0.27, K = 303.97##


#how did we do?
titpred<- titpop$pop[-length(titpop$pop)]*exp(0.79*(1-(titpop$pop[-length(titpop$pop)]/211)))
diff_titpred<- titpred-titpop$pop[-length(titpop$pop)]
diff_titpop <- titpop$pop[-1]-titpop$pop[-length(titpop$pop)]

plot(diff_titpop,diff_titpred)


########play around with missing data#############

titpopdf1<-as.data.frame(titpop$Broods)
colnames(titpopdf1) <- c('Broods')


###take away 15% of estimates##

missingdf1<-as.data.frame(lapply(titpopdf1, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.85, 0.15), size = length(cc), replace = TRUE) ]))

missingdf11<-missingdf1 %>% mutate(popplus1=Lag(titpop$Broods, -1))


#####
## NLS can handle missing data, but will accumulate bias##

Rickermodel22 <-nls(popplus1~Broods*exp(r*(1-(Broods/(1/a)))),start=list(r=0.1, a= 1/200), data=missingdf11)

### try with STAN###
missing_dat<- list(y=missingdf1$Broods, N=length(missingdf1$Broods))
missing_fit<- stan(file="ricker_.stan", data= missing_dat, iter = 1000, chains = 4)
