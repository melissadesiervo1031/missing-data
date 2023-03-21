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
pr <- dat %>% filter(site_name == id) %>% select(date, GPP, light, Q, GPP.upper,GPP.lower) %>% mutate(light.rel = light/max(light))
pr<-as.data.frame(pr)

pr1<-pr %>% select(GPP)


##quick plot of pr GPP dataset##
pineriverGPP <- ggplot(pr, aes(date, GPP))+
    geom_point(size = 2, color="chartreuse4") + 
    geom_line(size = 1, color="chartreuse4")+
    labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title="Pine River GPP")

##### #  ARIMA and STAN model with NO missing data #### ###note that there actually are some missing dates..ignore for now##

### ARIMA ###

X = matrix(c(pr$light.rel, pr$Q), ncol = 2)

arimafit <- Arima(pr$GPP, order = c(1,0,0), xreg = X)

#ar1:0.7139     intercept: 1.9943  light.rel= 5.2758  Q = -1.6564  #


### STAN #### 

# calculate observation error based on 95% CIs

sigma_obs <- (pr$GPP.upper - pr$GPP.lower)/3.92
pr <- pr %>%
  select(date, GPP, light, Q) %>%
  mutate(across(-date, ~zoo::na.approx(.)))%>%
  mutate(light.rel = light/max(light))


#Prep models ####
# model file: "model types/fixed_oi_light_centered.stan"
model_lq <- stan_model("GPP sim and real/Stan_code/AR1_light_Q_centered.stan")

#Create data object
data <- list(N = nrow(pr), P_obs = pr$GPP,
             sdo = sigma_obs, light = pr$light.rel, Q = pr$Q,
             miss_vec = rep(1, nrow(pr)))

#Run Stan
fit_lq <- rstan::sampling(object=model_lq, data = data,  
                          iter = 4000, chains = 4)

print(fit_lq,digits=5)

## beta1 (intercept) = 1.17345, beta2 = 2.93198 (light), beta 3 = -1.26293 (Q), phi = 0.62793, sdp = 1.17685 ###

## try again model w/ no error ###




####################Apply missing data function to Pine River  GPP ############

pine_missing_list<-makeMissing(timeSeries = pr$GPP, typeMissing = "random")


##add back in the date column and the covariates## 

pine_missing_list_11 <-lapply(X = pine_missing_list, FUN = function(X)   cbind.data.frame(GPP=X, date=as.Date(pr$date, format="%Y-%m-%d"), Q = pr$Q, light.rel=pr$light.rel))


###################### Data augmentaion in STAN #####################





################## Multiple imputations w/ AMELIA ###################

amelia1 <-lapply(X = pine_missing_list_11 , FUN = function(X)   amelia(X, ts="date", m=5, polytime=1))

### Amelia makes 5 version of each imputed dataset for each item in the list ###

##to look at one imputated dataset...##

m1try<-as.data.frame(amelia1[["propMissing_0.05"]][["imputations"]][["imp1"]])

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



