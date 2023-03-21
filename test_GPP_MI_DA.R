# Load packages
library(here)
library(tidyverse)
library(rstan)
library(Amelia)
library(forecast)
library(xts)
library(nlme)
library(tidyverse)

# Reference A. Stears' code with helpful function for removing data
# makeMissing()
source("missing_data_functions.R")

##upload pine river dataset##
dat <- read_csv('data/NWIS_MissingTS_subset.csv')
mdat <- read_csv('data/NWIS_MissingTSinfo_subset.csv')

id <- mdat$site_name[4]
pr <- dat %>% filter(site_name == id) %>% select(date, GPP, light, Q) %>% mutate(light.rel = light/max(light))
pr1<-pr %>% select(GPP)


##quick plot of pr GPP dataset##
pineriverGPP <- ggplot(pr, aes(date, GPP))+
    geom_point(size = 2, color="chartreuse4") + 
    geom_line(size = 1, color="chartreuse4")+
    labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title="Pine River GPP")


##### make 10 copies of the pine river dataset into a list##

pine_listwdate<-replicate(n = 10, expr = {data.frame(pr)}, simplify = F)

pine_list<-replicate(n = 10, expr = {data.frame(pr1)}, simplify = F)

###

####################Apply missing data function to list of PR GPP ############

pine_missing_list <-lapply(X = pine_listwdate, FUN = function(X)   makeMissing(timeSeries = X$GPP, typeMissing = "random")  )

pine_missing_list_1<-pine_missing_list[[1]]  


## data augmention ##



################## Multiple imputations w/ AMELIA ###################

##need to add back in the date column and the covariates## 

pine_missing_list_11 <-lapply(X = pine_missing_list_1, FUN = function(X)   cbind(X, day=seq(1, 360), Q = pr$Q, light.rel=pr$light.rel))

amelia1 <-lapply(X = pine_missing_list_11, FUN = function(X)   amelia(X, ts="day", m=5, polytime=1))

### Amelia makes 5 version of each imputed dataset for each item in the list ###

##to look at one imputated dataset...##

m1try<-as.data.frame(amelia1[["propMissing_0.05"]][["imputations"]][["imp1"]])

m1try

##to look at 1 group of 5 imputed dataset ## 

m1to5<-amelia1[["propMissing_0.05"]][["imputations"]]


### ARIMA model to run on 1 imputed datasets###

X = matrix(c(pr$light.rel, pr$Q), ncol = 2)

fit <- Arima(xts(m1try$X, order.by = as.Date(m1try$day)), order = c(1,0,0), xreg = X)

###
modelOutput_list <- replicate(length(5), rep(NULL), simplify = FALSE)

### not yet working but close...##

for(i in 1:5){
  ## fit an ARIMA model to each AMELIA imputed dataset ##
  arimamissing_i <- Arima(xts(m1to5[[i]][,1], order.by = as.Date(m1to5[[i]][,2])), order = c(1,0,0), xreg = X)
  modelOutput_list[[i]] <- summary(arimamissing_i)$coefficients
}



