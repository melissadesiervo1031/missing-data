
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(plyr)
library(shinystan)
library(faux)
library(bayesplot)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(data.table)
library(Amelia)

## Simulate data
N<-365 #length of data
t<-1:N #time
time<-seq(as.POSIXct("2020/1/1"),as.POSIXct("2020/12/30"), by="day") #for light

##GPP based on time-series model with known paramters
set.seed(553)
x<-NA
sdp <- 0.1
phi<- 0.5
x[1]<-0.1
b1<-0.5





##light
# From Yard et al. (1995) Ecological Modelling.  Remember your trig?  
# calculate light as umol photon m-2 s-1.
# Arguments are:  
# time = a date and time input (posixct object)
# lat = latitude of field site
# longobs = longitude of field site
# longstd = standard longitude of the field site (NOTE: watch daylight savings time!!!). For PST, longstd is be 120 degrees. But during PDT it is 105 degrees. MST is 105 deg. MDT is 90. 


# convert degrees to radians
radi<-function(degrees){(degrees*pi/180)}

# function to estimate light
lightest<- function (time, lat, longobs, longstd) {
  jday<-yday(time)
  E<- 9.87*sin(radi((720*(jday-81))/365)) - 7.53*cos(radi((360*(jday-81))/365)) - 1.5*sin(radi((360*(jday-81))/365))
  LST<-as.numeric(time-trunc(time))
  ST<-LST+(3.989/1440)*(longstd-longobs)+E/1440
  solardel<- 23.439*sin(radi(360*((283+jday)/365)))
  hourangle<-(0.5-ST)*360
  theta<- acos( sin(radi(solardel)) * sin(radi(lat)) + cos(radi(solardel)) * cos(radi(lat)) * cos(radi(hourangle)) )
  suncos<-ifelse(cos(theta)<0, 0, cos(theta))
  GI<- suncos*2326
  GI	
  
}
light<-lightest(time, 47.8762, -114.03, 105) #Flathead Bio Station just for fun
light.l<-log(light)
light.c<-(light-mean(light))/sd(light)
light.rel<-light/max(light)
mean(light.c)
sd(light.c)



## Set the seed, so we can reproduce the results
set.seed(553)
## For-loop that simulates the state through time, using i instead of t,
for (t in 2:N){
  x[t] = x[t-1]*phi+b1*light.rel[t]+rnorm(1, 0, sdp)
}

## plot simulated data
plot(1:N, x[1:N],
     pch=3, cex = 0.8, col="blue", ty="o", lty = 3,
     xlab = "t", ylab = expression(x[t]),
     xlim = c(0,N), ylim = c(min(x), max(x)),
     las = 1)


##Loads lists of missing and observed data lengths from 10 yrs of data collected at Shatto Ditch
obs<-readRDS("C:/Users/matt/IDrive-Sync/Postdoc/Estimating missing data/daily_predictions/SDW_obs_length.RDS")
miss<-readRDS("C:/Users/matt/IDrive-Sync/Postdoc/Estimating missing data/daily_predictions/SDW_miss_length.RDS")


##Create many simulated datasets with varying amounts of missing data 
n.datasets<-10
y_miss<-list(x)
y_miss<-rep(y_miss[1], times=n.datasets)


##set up loop to create missing datasets from y_miss. 
#The loop pulls random lengths of observed and missing data and makes missing data zero. 
#It switches between pulling a missing length and replaces those values with zero,
#then it picks an observed length and jumps ahead (keeping the data points in the jump).
#The loop for each individual dataset ends when it exceeds the length N
#The end result are N=n.datasets with varying missing data gaps and percents of missing data

#Create objects for the loop
sample.miss<-list(NA)
sample.obs<-list(NA)
sample.miss<-rep(sample.miss, times=n.datasets)
sample.obs<-rep(sample.obs, times=n.datasets)
set.seed(500)

#For loop for N=n.datasets. The end result are N=n.datasets with 0 for missing data
for (i in 2:n.datasets){
  #samples 200 times from SDW lengths with replacement
  sample.miss[[i]]<-sample(miss,200,replace=TRUE)
  sample.obs[[i]]<-sample(obs,200,replace=TRUE)
  
  #Initiates new dataset. It will alwasys start with some length of observed data. 
  new<-sample.obs[[i]][1]
  
  #For loop for each individual dataset
  for(j in 1:365){
    
    remove<-seq(from=new+1, to=new+sample.miss[[i]][j], by=1)
    y_miss[[i]][remove]<-0
    new<-new+sample.miss[[i]][j]+sample.obs[[i]][j+1]
    if (new>=N){
      break
    }
  }
}

#Assign NAs to the missing data in each list
for(i in 1:n.datasets){
  z<-y_miss[[i]]
  z[which(z==0)]<-NA
  y_miss[[i]]<-z
}

#Cuts each down to size, as they will often end with missing data
for (i in 1:n.datasets){
  y_miss[[i]]<-y_miss[[i]][1:365]
}


#Check to make sure
summary(y_miss)


#Calculate the average length of data gaps in each dataset
## Calculate length of data gaps
# Prepare vectors for loop
miss_vec<-list(NA)
miss_vec<-rep(miss_vec, times=n.datasets)
length_vec<-list(NA)
length_vec<-rep(length_vec, times=n.datasets)

# Create dummy variable for present or not
dummy<-list(NA)
dummy<-rep(dummy, times=n.datasets)

for (i in 1:n.datasets){
  dummy[[i]]<-ifelse(is.na(y_miss[[i]]),1,0)
}

#loop over data and return the length of missing data gaps
for (i in 1:n.datasets){
  # First data point will always be observed
  miss_vec[[i]][1]<-dummy[[i]][1]
  length_vec[[i]][1]<-1
  
  for (j in 1:N){
    miss_vec[[i]][j]<-ifelse(dummy[[i]][j]==1,1+miss_vec[[i]][j-1], 0) #If data is missing add to yesterdays value by 1. If data observed then 0.
    length_vec[[i]][[j]]<-ifelse(miss_vec[[i]][j]==0,miss_vec[[i]][j-1], NA) #Save length of gap on the day data is observed (0) after a gap
    # length_miss[i]<-ifelse(length_miss[i]==0,NA,length_miss[i]) #Replace days with observed data with NA
  }
}

##Calculate mean of missing data gaps
missing<-lapply(length_vec,na.omit)
missing<-lapply(missing, function(missing) {missing[missing!=0]})
mean.missing.gap<-unlist(lapply(missing,mean))
mean.missing.gap

##first dataset has no gaps (is NaN)
mean.missing.gap[[1]]<-0
hist(mean.missing.gap)
#Calculate the percent of missing data
miss.num<-sapply(y_miss, function(x) sum(length(which(is.na(x)))))
prop.miss<-round(miss.num/length(x)*100)
prop.miss
hist(prop.miss)

##Create objects with the location (index) of missing or observed data AND number
## of missing or observed data
#FOR SINGLE VECTOR
#y_index_mis <- which(is.na(y_miss)) # Identify the rows for each variable in airquality with NA
#y_index_obs <- which(!is.na(y_miss)) # Identify the rows for each variable in airquality with observed data
#y_nMiss <- length(y_index_mis)# How many NAs?
#y_nObs <- length(y_index_obs) # How many NAs?

##Create objects with the location (index) of missing or observed data AND number
## of missing or observed data
#FOR LIST OF VECTORS
y_index_mis <- lapply(y_miss,function(var){which(is.na(var))}) # Identify the rows for each variable in airquality with NA
y_index_obs <- lapply(y_miss,function(var){which(!is.na(var))}) # Identify the rows for each variable in airquality with observed data
y_nMiss <- lapply(y_index_mis,function(var){length(var)}) # How many NAs?
y_nObs <- lapply(y_index_obs,function(var){length(var)}) # How many NAs?


##Replace NAs with arbitrary number to make Stan happy
for(i in 1:n.datasets){
  r<-y_miss[[i]]
  r[y_index_mis[[i]]]<--100 #arbitrary number
  y_miss[[i]]<-r
} 

summary(y_miss)



sink("mcelreath_miss_ss.stan")

cat("
 

/*----------------------- Functions --------------------------*/  
  functions{
    vector merge_missing( int[] y_index_mis, vector y_miss, vector y_imp) {
    int N = dims(y_miss)[1];
    int N_miss = dims(y_imp)[1];
    vector[N] merged;
    merged = y_miss;
    for (i in 1:N_miss)
        merged[y_index_mis[i] ] =y_imp[i];
    return merged;    
    }
  }
  
/*----------------------- Data --------------------------*/
  /* Data block: defines the objects that will be inputted as data */
  data {
    int N; // Length of state and observation time series
    vector[N] y_miss; // Observations
    real z0; // Initial state value
    int y_nMiss; // number of missing values
    int y_index_mis[y_nMiss]; // index or location of missing values within the dataset
  }
  
/*----------------------- Parameters --------------------------*/
  /* Parameter block: defines the variables that will be sampled */
  parameters {
    vector<lower = 0>[y_nMiss] y_imp;// Missing data
    real<lower=0> sdp; // Standard deviation of the process equation
    //real<lower=0> sdo; // Standard deviation of the observation equation
    real b0;
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
   // vector[N] z; // State time series
  }
  

  /*----------------------- Model --------------------------*/
  /* Model block: defines the model */
  model {
    //merge imputed and observed data
    vector[N] B_merge;
    B_merge= merge_missing(y_index_mis, to_vector(y_miss), y_imp);
   
    // Prior distributions
    //sdo ~ normal(0, 1);
    sdp ~ normal(0, 1);
    phi ~ beta(1,1);
    b0 ~ normal(0,5);
  
    
    // Distribution for the first state
  B_merge[1] ~ normal(z0, sdp);
   //y_miss[1]~ normal(B_merge[1],sdo);
   
    // Distributions for all other states
    for(t in 2:N){
       B_merge[t] ~ normal(b0+ B_merge[t-1]*phi, sdp);
       //y_miss[t] ~ normal(B_merge[t], sdo); // observation model with fixed observation error
    
    }
   
  }
  
  
    "
    ,fill=TRUE)
sink()
closeAllConnections()


fit.miss <- vector("list",n.datasets)
model<-"mcelreath_miss_ss.stan"
model<-stan_model(model)
for(i in 1:10){
  ##Load data
  data <- list(   N = length(y_miss[[i]]),
  y_nMiss = unlist(y_nMiss[[i]]),
  y_index_mis = unlist(y_index_mis[[i]]),
  y_miss= unlist(y_miss[[i]]),
  z0=y_miss[[i]][1])
  fit.miss[[i]]<- rstan::sampling(object=model, data = data,  iter = 4000, chains = 4, control=list(max_treedepth = 15))
  print(i)
  rstan::check_hmc_diagnostics(fit.miss[[i]])
}

traceplot(fit.miss[[1]], pars=c("phi","sdp","b0"))


##Pull param estimates into list
fit_summary_pars_bayes <- vector("list",10)
for (i in 1:10){
  fit_summary_pars_bayes[[i]]<-(summary(fit.miss[[i]], pars=c("sdp","phi", "b0"), probs=c(0.025,.5,.975))$summary)
}
fit_summary_pars_bayes[[2]]

##Unlist,cleanup, and add factors
fit_summary_bayes<-as.data.frame(do.call(rbind.data.frame, fit_summary_pars_bayes)) #Unlist
fit_summary_bayes$param<-rep(c("sdp","phi","b0"), times=10 )#add parameter factor
fit_summary_bayes$prop.missing<-rep(prop.miss[1:10], times=rep(3, times=10)) #add prop of missing data
fit_summary_bayes$gap.missing<-rep(mean.missing.gap[1:10], times=rep(3, times=10)) #add prop of missing data
row.names(fit_summary_bayes)<-1:(10*3) #remove row names
fit_summary_bayes$param<-as.factor(fit_summary_bayes$param) #make factor
#fit_summary_bayes$prop.missing<-as.factor(fit_summary_bayes$prop.missing) #make factor
colnames(fit_summary_bayes)<-c("mean", "se-mean", "sd", "min", "med", "high","n_eff", "rhat", "param","prop.missing", "gap.missing")
summary(fit_summary_bayes) #check it looks good
head(fit_summary_bayes)

fit_summary_bayes<-fit_summary_bayes[order(fit_summary_bayes$prop.missing),]
known.data<-c(sdp,phi,0)
known<-as.data.frame(known.data)

known.param<-c("sdp","phi","b0")
#known.missing<-c(5, 5, 5,5)
known<-as.data.frame(cbind(known.data, known.param))
known$known.data<-as.numeric(as.character(known$known.data))

ggplot(data=fit_summary_bayes, aes(x=mean, y=prop.missing ))+
  geom_point(aes( color=param, group=param),size=3)+
  theme_classic()+
  theme(legend.position="top")+
  ylab("Percent of Missing Data")+
  xlab("Parameter Estimate")+
  scale_color_manual(values=c("blue", "green", "black"))+
  #scale_x_continuous(limits=c(0,1))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  geom_vline(xintercept = known.data,color=c("black", "green", "blue"))



##Calculate denominator of rwci for easier calculation
fit_summary_bayes$rwci.den<-rep(fit_summary_bayes$high[1:3]-fit_summary_bayes$min[1:3], times=10)
fit_summary_bayes$known.param<-rep(c(0.1, 0.5, 0), times=10)

##Calculate bias (difference from known) and rwci (relative width of the credible interval)

diff<-ddply(fit_summary_bayes, c("prop.missing", "param"),summarize,
            bias=mean-known.param, 
            rwci= (high-min)/rwci.den
)

fit_summary_bayes$bias<-diff$bias
fit_summary_bayes$rwci<-diff$rwci
known.param<-as.data.frame(cbind(as.numeric(fit_summary_bayes$known.param[1:3]), as.character(fit_summary_bayes$param[1:3])))
colnames(known.param)<-c("known","param")
known.param$known<-as.numeric(known.param$known)

ggplot(data=fit_summary_bayes, aes(x=mean, y=prop.missing ))+
  geom_point(aes( color=param, group=param),size=3)+
  theme_classic()+
  theme(legend.position="right")+
  xlab("Parameter value")+
  ylab("Percent of Missing Data")+
  scale_color_manual(values=c("blue", "green", "black"))+
  #scale_x_continuous(limits=c(-1,1))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  geom_vline(xintercept =known.param$known,color=c("black", "green", "blue"))

ggplot(data=fit_summary_bayes, aes(x=log(rwci), y=prop.missing ))+
  geom_point(aes( color=param, group=param),size=3)+
  theme_classic()+
  theme(legend.position="right")+
  ylab("Percent of Missing Data")+
  xlab("log Relative width of 95% credible interval log")+
  scale_color_manual(values=c("blue", "green", "darkgray", "black"))+
  #scale_x_continuous(limits=c(-.05,1.1))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  geom_vline(xintercept = 0,color=c("black"))

ggplot(data=fit_summary_bayes, aes(x=bias, y=gap.missing ))+
  geom_point(aes( color=param, group=param),size=3)+
  theme_classic()+
  theme(legend.position="right")+
  xlab("Parameter value")+
  ylab("Percent of Missing Data")+
  scale_color_manual(values=c("blue", "green", "darkgray", "black"))+
  #scale_x_continuous(limits=c(-1,1))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  geom_vline(xintercept =known.param$known,color=c("black","darkgray", "blue", "green" ))


