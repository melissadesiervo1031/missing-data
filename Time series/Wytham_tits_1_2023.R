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


###calculate lambda between years###

##add variable t+1##

titpop$poptplus1 <- Lag(titpop$Broods, -1)

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

#Fairly close match with STAN and NLS### ## r =  0.491790068 , alpha = 0.001561374 (K = 315)




#### solve for estimates with ARIMA ###

##helpful links https://ionides.github.io/531w16/final_project/Project19/final.html ##
#https://towardsdatascience.com/state-space-model-and-kalman-filter-for-time-series-prediction-basic-structural-dynamic-linear-2421d7b49fa6#
#https://lbelzile.github.io/timeseRies/state-space-models-and-the-kalman-filter.html



########play around with missing data#############

titpopdf1<-as.data.frame(titpop$Broods)
colnames(titpopdf1) <- c('Broods')


###create vector of % of missing data ###

n_missing<-seq(from = 0, to = 70, by = 5)


## duplicate tit pop dataset in a list to use forloop to create missing ###

listmissingtitdf <- replicate(length(n_missing),titpopdf1, simplify = FALSE)


###duplicate tit pop into rows of a dataframe with missing values###

missingtitdf<-data.frame(replicate(length(n_missing),titpopdf1))

#for loop to make TITPOP dataset with varying levels of missingness###


##not quite right###
for (i in 1:length(listmissingtitdf)) {
  
  listmissingtitdf[[i]]<-sample(c(TRUE, NA), prob = c(0.95, 0.05), size = length(titpopdf1$Broods), replace = TRUE)
  
}


##not quite right###
for (i in 1:length(listmissingtitdf)) {
  
  listmissingtitdf[[i]]<-sample(c(TRUE, NA), prob = c(0.95, 0.05), size = length(titpopdf1$Broods), replace = TRUE)
  
}

##removes 20% of data from each column###
missingdftry<-missingtitdf  %>% mutate_all(~ifelse(sample(c(TRUE, FALSE), size = length(.), replace = TRUE, prob = c(0.8, 0.2)),
                            as.character(.), NA))


for (i in 1:length(listmissingtitdf)) {
  
 yesno[[i]]<-sample(c(1, NA), prob = c(0.95, 0.05), size = length(titpopdf1$Broods), replace = TRUE)
 listmissingtitdf[[i]]<-listmissingtitdf[[i]]*yesno[[i]]
   
}


###take away 5% of estimates##

missingdf5<-as.data.frame(lapply(titpopdf1, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.95, 0.05), size = length(cc), replace = TRUE) ]))
missingdf55<-missingdf5 %>% mutate(popplus1=Lag(missingdf5$Broods, -1))


###take away 20% of estimates##

missingdf20<-as.data.frame(lapply(titpopdf1, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.80, 0.20), size = length(cc), replace = TRUE) ]))
missingdf20<-missingdf20 %>% mutate(popplus1=Lag(missingdf20$Broods, -1))


###take away 40% of estimates##

missingdf40<-as.data.frame(lapply(titpopdf1, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.60, 0.40), size = length(cc), replace = TRUE) ]))
missingdf40<-missingdf40 %>% mutate(popplus1=Lag(missingdf40$Broods, -1))


###take away 70% of estimates##

missingdf70<-as.data.frame(lapply(titpopdf1, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.30, 0.70), size = length(cc), replace = TRUE) ]))
missingdf70<-missingdf70%>% mutate(popplus1=Lag(missingdf70$Broods, -1))



###take away 90% of estimates##

missingdf90<-as.data.frame(lapply(titpopdf1, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.10, 0.90), size = length(cc), replace = TRUE) ]))
missingdf90<-missingdf90%>% mutate(popplus1=Lag(missingdf90$Broods, -1))




#####
## NLS can handle missing data (it removes them) but it will accumulate bias##

###""TRUE""" estimates from no missing data### 

Rickermodelfull <-nls(poptplus1~Broods*exp(r-alpha*Broods),start=list(r=0.1, alpha= 1/200), data=titpop)

## r = 0.3083 , K = 308## (alpha = 0.000997)

Rickermissing5 <-nls(popplus1~Broods*exp(r-alpha*Broods),start=list(r=0.1, alpha= 1/200), data=missingdf55)
#na.action(Rickermodel3) ## tells us what rows were omitted ###

Rickermissing20 <-nls(popplus1~Broods*exp(r-alpha*Broods),start=list(r=0.1, alpha= 1/200), data=missingdf20)

Rickermissing40 <-nls(popplus1~Broods*exp(r-alpha*Broods),start=list(r=0.1, alpha= 1/200), data=missingdf40)

Rickermissing70 <-nls(popplus1~Broods*exp(r-alpha*Broods),start=list(r=0.1, alpha= 1/200), data=missingdf70)


##doesn't work###

Rickermissing90 <-nls(popplus1~Broods*exp(r-alpha*Broods),start=list(r=0.1, alpha= 1/200), data=missingdf90)

### pull out estimates ###

Rickermodelfullest<-as.data.frame(summary(Rickermodelfull )$coefficients) 

Rickermissingfullestdf<-Rickermodelfullest %>% mutate(missing=0) %>% tibble::rownames_to_column("param")


Rickermissing5est<-as.data.frame(summary(Rickermissing5)$coefficients) 

Rickermissing5estdf<-Rickermissing5est %>% mutate(missing=5) %>% tibble::rownames_to_column("param")


Rickermissing20est<-as.data.frame(summary(Rickermissing20)$coefficients) 

Rickermissing20estdf<-Rickermissing20est %>% mutate(missing=20) %>% tibble::rownames_to_column("param")


Rickermissing40est<-as.data.frame(summary(Rickermissing40)$coefficients) 

Rickermissing40estdf<-Rickermissing40est %>% mutate(missing=40) %>% tibble::rownames_to_column("param")


Rickermissing70est<-as.data.frame(summary(Rickermissing70)$coefficients) 

Rickermissing70estdf<-Rickermissing70est %>% mutate(missing=70) %>% tibble::rownames_to_column("param")


##merge the dataframe w/ estimates from models with 0 - 40 % missing data ##

Rickermissingalldf<-rbind(Rickermissingfullestdf, Rickermissing5estdf, Rickermissing20estdf, Rickermissing40estdf, Rickermissing70estdf)



###plot Std error of estimates over missing-ness##

ggplot(data=Rickermissingalldf, aes(x=missing, y=Estimate))+
  facet_wrap(~param, scales="free")+
  geom_errorbar(aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`), width=.2)+
  geom_point(size=3)+
  theme_classic()+
  xlab("Percent of Missing Data")+
  ylab("Parameter Estimate")
  

