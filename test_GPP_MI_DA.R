# Load packages
library(here)
library(tidyverse)
library(rstan)
library(Amelia)

# Reference A. Stears' code with helpful function for removing data
# makeMissing()
source("missing_data_functions.R")

##upload pine river dataset##
dat <- read_csv('data/NWIS_MissingTS_subset.csv')
mdat <- read_csv('data/NWIS_MissingTSinfo_subset.csv')

id <- mdat$site_name[4]
pr <- dat %>% filter(site_name == id) %>% select(date, GPP)
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

################## Multiple imputations w/ AMELIA ###################

##need to add back in the date column## 

pine_missing_list_11 <-lapply(X = pine_missing_list_1, FUN = function(X)   cbind(X, day=seq(1, 360)))

amelia1 <-lapply(X = pine_missing_list_11, FUN = function(X)   amelia(X, ts="day", m=5, polytime=1))

### Amelia makes 5 version of each imputed dataset for each item in the list ###

##to look at one imputation...##

m1try<-amelia1[["propMissing_0.05"]][["imputations"]][["imp1"]]



