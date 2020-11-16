##Load libraries

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(dplyr)
library(shinystan)
library(tidyverse)
library(reshape2)

##Load airquality dataset

data(airquality, package = "datasets")
temp<-airquality$Temp

##Make date column

airquality$date<- paste(airquality$Month,airquality$Day, sep="." )
airquality$date<-as.Date(airquality$date, "%m.%d")

## Plot data
ggplot(data=airquality, aes(x=date, y=Temp))+
  geom_point( size=2)+
  theme_classic()+
  theme(legend.position="top")+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))

##Make a copy of dataset
airquality2<-airquality

##Randomly remove temperature data and replace with NA

set.seed(568)
temp_na<- which(temp %in% sample(temp,5)) ##2,5,15
airquality2$Temp[temp_na]<-NA

##Randomly remove 3 continuous 5-day chunks
#set.seed(568)
#temp_na<-(which(temp %in% sample(temp, 1.5)))
#n=5




##Create objects with the location index of missing and observed data AND 
#number of missing and observed data

missing <- lapply(airquality2,function(var){which(is.na(var))}) # Identify the rows for each variable in airquality with NA
observed <- lapply(airquality2,function(var){which(!is.na(var))}) # Identify the rows for each variable in airquality with observed data
nMissing <- lapply(missing,function(var){length(var)}) # How many NAs?
nObserved <- lapply(observed,function(var){length(var)}) # How many NAs?

##Rename column headers for clarity

names(missing) <- paste0(names(missing),'_index_mis') # Rename indices to variable_index_mis
names(nMissing) <- paste0(names(nMissing),'_nMiss') # Rename n_miss to variable_nMiss
names(observed) <- paste0(names(observed),'_index_obs') # Rename indices to variable_index_obs
names(nObserved) <- paste0(names(nObserved),'_nObs') # Rename n_Obs to variable_nObs

##Replace NAs with arbitrary number to make Stan happy

airquality2[is.na(airquality2)] <- -100

##Stan Code

sink("temp_AR_blog_obs.stan")

cat("
    data {
    int N; #Number of observations
    vector[N] temp; #Response variable, including missing values
    int Temp_nMiss;
    int Temp_index_mis[Temp_nMiss];
    }
    
    parameters {
    vector[Temp_nMiss] temp_imp;//Missing data
    real<lower = 0> sigma; //  standard deviation
     real alpha; // Constant
    real beta;  // AR
    }
    
    transformed parameters { 
    vector[N] y;
    y = temp; 
    y[Temp_index_mis] =temp_imp;
    } 
    
    model {
    for (n in 2:N){
    y[n] ~ normal(alpha+beta*y[n-1], sigma);
      }
    sigma ~ normal(1, 1);
    }
    
    "
    ,fill=TRUE)
sink()
closeAllConnections()

##Load data

aq.temp <- list(N = length(temp),
                Temp_nMiss = nMissing$Temp_nMiss,
                Temp_nObs = nObserved$Temp_nObs,
                Temp_index_mis = missing$Temp_index_mis,
                Temp_index_obs = observed$Temp_index_obs,
                temp = airquality2$Temp
                )

##Run Stan

fit <- stan("temp_AR_blog_obs.stan", data = aq.temp,  iter = 1000, chains = 4)

##basic fit evaluation and extraction

print(fit)
class(fit)
traceplot(fit, pars="sigma")
fit_extract<-rstan::extract(fit)

##Create object with estimated missing data

fit_summary<-summary(fit, probs=c(0.05,.5,.95))$summary
my_data <- as_tibble(fit_summary)
y_est_data<-my_data %>% slice(1:length(missing$Temp_index_mis))

##Create object with observed data

y_miss_index<-missing$Temp_index_mis
y_obs_data<-as.data.frame(cbind(airquality$date[y_miss_index],airquality$Temp[y_miss_index]))
colnames(y_obs_data)<-c("date","temp")
y_obs_data$date<-as.Date(y_obs_data$date, origin = "1969-12-30")

##Merge observed and estimated

y_combined<-cbind(y_obs_data,y_est_data)

##rename column names for clarity

colnames(y_combined)<-c("date", "observed", "estimated", "se_mean_est", "sd_est", "low", "median","high", "n_eff", "rhat")

##check that column names are correct

head(y_combined)

##Plot observed and estimated direct comparison

ggplot(data=y_combined, aes(x=observed, y=estimated))+
  geom_point( size=3)+
  geom_abline(intercept=0, slope=1)+
  theme_classic()+
  geom_errorbar(aes(ymin=low, ymax=high), width=0.5)+
  theme(legend.position="top")+
  ylab("Estimated Temperature")+
  xlab("Observed Temperature")+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))

##Plot observed and estimated time series

ggplot(data=y_combined, aes(x=date, y=estimated))+
  geom_line(data=airquality, aes(x=date, y=Temp, color="black"), size=1)+
  geom_point(aes(color="red"), size=3)+
  geom_errorbar(data=y_combined,aes(x=date,ymin=low, ymax=high))+ 
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  scale_color_identity(guide = "legend", name=element_blank(), labels=c("Raw data", "Estimated points"))+
  theme(legend.position="top")+
  ylab("Temperature")
