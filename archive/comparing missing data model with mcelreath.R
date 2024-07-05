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
library(rethinking)
###Simulate data
N<-365 #length of data
t<-1:N #time
time<-seq(as.POSIXct("2020/1/1"),as.POSIXct("2020/12/30"), by="day") #for light

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
mean(light.c)
sd(light.c)


##GPP based on time-series model with known paramters
set.seed(553)
x<-NA
sdp <- 0.01
phi <-0.8
b0<-0.1
b1<-0.1
x[1]<-3.5
# Set the seed, so we can reproduce the results
set.seed(553)
# For-loop that simulates the state through time, using i instead of t,
for (t in 2:N){
  x[t] = b0+phi*x[t-1]+light.l[t]*b1+rnorm(1, 0, sdp)
}

#Estimate latent state
#for (t in 2:N){
#      x[t] = b0+phi*x[t-1]+rnorm(1, 0, sd.p)
#}


#Estimate observed data
#set.seed(553)
#sd.o<-0.2
#sdo <-rnorm(N, 0, sd.o)
#y <-x + sdo

#sdo_sd<-mean(abs(sdo))

y_full<-as.data.frame(cbind(x,time))

ggplot(data=y_full, aes(x=time, y=x))+
  geom_point(color="black", size=3)+
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  theme(legend.position="top")+
  ylab("Data")+
  xlab("Time (days)")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle("Full dataset")

t<-1:N
y_miss_40<-as.data.frame(cbind(y_miss[[6]],t))

ggplot(data=y_miss_40, aes(x=t, y=V1))+
  geom_point(color="black", size=3)+
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  theme(legend.position="top")+
  ylab("Data")+
  xlab("Time (days)")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle("60% removed")



##Plot observed and latent
plot(1:N, x[1:N],
     pch=3, cex = 0.8, col="blue", ty="o", lty = 3,
     xlab = "t", ylab = expression(z[t]),
     xlim = c(0,N), ylim = c(min(x), max(x)),
     las = 1)
points(1:N, y,
       pch = 19, cex = 0.7, col="red", ty = "o")
legend("top",
       legend = c("Obs.", "Latent states"),
       pch = c(3, 19),
       col = c("blue", "red"),
       lty = c(3, 1),
       horiz=TRUE, bty="n", cex=0.9)


##Force some missing data-BY DAY##
#Note-the model will still work if no data is missing but the objects in the next block of code (i.e., "y_index_mis", etc.)
#still need to be created

set.seed(50)

##Load data from SDW characterizing lengths
obs<-readRDS("C:/Users/matt/IDrive-Sync/Postdoc/Estimating missing data/daily_predictions/SDW_obs_length.RDS")
miss<-readRDS("C:/Users/matt/IDrive-Sync/Postdoc/Estimating missing data/daily_predictions/SDW_miss_length.RDS")

##Replicate simulated data
n.datasets<-5
y_miss<-list(x)
y_miss<-rep(y_miss[1], times=n.datasets)


##set up loop to create missing datasets from y_miss. 
##The loop pulls random lengths of observed and missing data and makes missing data zero. 
##It switches between pulling a missing length and replaces those values with zero,
##then it picks an observed length and jumps ahead (keeping the data points in the jump).
##The loop for each individual dataset ends when it exceeds the length N
##The end result are N=n.datasets with varying missing data gaps and percents of missing data

#Create objects for the loop
sample.miss<-list(NA)
sample.obs<-list(NA)
sample.miss<-rep(sample.miss, times=n.datasets)
sample.obs<-rep(sample.obs, times=n.datasets)

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


##Stan Code
sink("mcelreath_miss.stan")

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
    vector[N] light;
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
    real b1;
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    //vector[N] z; // State time series
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
    b1 ~ normal(0,5);
    
    // Distribution for the first state
   B_merge[1] ~ normal(z0, sdp);
    
    // Distributions for all other states
    for(t in 2:N){
       B_merge[t] ~ normal(b0+ B_merge [t-1]*phi+light[t]*b1, sdp);
    }
    // Distributions for the observations
   //for(t in 1:N){
    
    //  y[t] ~ normal(z[t], sdo);
   //}
  }
  
  
    "
    ,fill=TRUE)
sink()
closeAllConnections()

fit.mcelreath.miss <- vector("list",n.datasets)
model<-"mcelreath_miss.stan"
model<-stan_model(model)
for(i in 1:5){
  ##Load data
  data <- list(   N = length(y_miss[[i]]),
                  y_nMiss = unlist(y_nMiss[[i]]),
                  y_index_mis = unlist(y_index_mis[[i]]),
                  y_miss= unlist(y_miss[[i]]),
                  light=light.l,
                  z0=y_miss[[i]][1])
  
  ##Run Stan
  fit.mcelreath.miss[[i]]<- rstan::sampling(object=model, data = data,  iter = 4000, chains = 4, control=list(max_treedepth = 15))
  print(i)
}

fit_summary_pars_mcelreath<-as.data.frame((summary(fit.mcelreath.miss[[2]], pars=c("sdp","phi", "b1", "b0"), probs=c(0.025,.5,.975))$summary))

fit_summary_pars_mcelreath <- vector("list",5)
for (i in 1:5){
  fit_summary_pars_mcelreath[[i]]<-(summary(fit.mcelreath.miss[[i]], pars=c("sdp","phi", "b1", "b0"), probs=c(0.025,.5,.975))$summary)
}
fit_summary_pars_mcelreath[[5]]
#### model 2= ~340 seconds





##Stan Code
sink("matt.miss.stan")

cat("
 
/*----------------------- Data --------------------------*/
  /* Data block: defines the objects that will be inputted as data */
  data {
    int N; // Length of state and observation time series
    vector[N] y_miss; // Observations
    real z0; // Initial state value
    vector[N] light;
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
    real b1;
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    //vector[N] z; // State time series
  }
  transformed parameters { 
    vector[N] y;
    y=y_miss; // makes the data a transformed variable
    y[y_index_mis] =y_imp; // replaces missing data in y with estimated data
    } 
  /*----------------------- Model --------------------------*/
  /* Model block: defines the model */
  model {
    // Prior distributions
    //sdo ~ normal(0, 1);
    sdp ~ normal(0, 1);
    phi ~ beta(1,1);
    b0 ~ normal(0,5);
    b1 ~ normal(0,5);
    // Distribution for the first state
    y[1] ~ normal(z0, sdp);
    // Distributions for all other states
    for(t in 2:N){
      y[t] ~ normal(b0+y[t-1]*phi+light[t]*b1, sdp);
    }
    // Distributions for the observations
   //for(t in 1:N){
    
    //  y[t] ~ normal(z[t], sdo);
   //}
  }
  generated quantities {
 vector[N] y_rep; // replications from posterior predictive dist
 y_rep[1]=normal_rng(y[1], 0.1);
 
 for (t in 2:N) {
 y_rep[t]=normal_rng(b0+y[t-1]*phi+light[t]*b1, sdp);
 }
 
  }

    "
    ,fill=TRUE)
sink()
closeAllConnections()



fit.matt.miss <- vector("list",n.datasets)
model<-"matt.miss.stan"
model<-stan_model(model)
for(i in 1:5){
  ##Load data
  data <- list(   N = length(y_miss[[i]]),
                  y_nMiss = unlist(y_nMiss[[i]]),
                  y_index_mis = unlist(y_index_mis[[i]]),
                  y_miss= unlist(y_miss[[i]]),
                  light=light.l,
                  z0=y_miss[[i]][1])
  
  ##Run Stan
  fit.matt.miss[[i]]<- rstan::sampling(object=model, data = data,  iter = 4000, chains = 4, control=list(max_treedepth = 15))
  print(i)
}
##model 2= ~ 360 seconds 


fit_summary_pars_matt<-as.data.frame((summary(fit.matt.miss[[2]], pars=c("sdp","phi", "b1", "b0"), probs=c(0.025,.5,.975))$summary))

fit_summary_pars_matt <- vector("list",5)
for (i in 1:5){
  fit_summary_pars_matt[[i]]<-(summary(fit.matt.miss[[i]], pars=c("sdp","phi", "b1", "b0"), probs=c(0.025,.5,.975))$summary)
}


for (i in 1:5){
plot(as.data.frame(fit_summary_pars_matt[[i]])$mean~as.data.frame(fit_summary_pars_mcelreath[[i]])$mean)
  abline(a=0, b=1)
}
abline(a=0, b=1)


fit_summary_matt<-summary(fit.matt.miss[[2]], probs=c(0.025,.5,.975))$summary

my_data <- as_tibble(fit_summary_matt)
y_est_data_matt<-my_data %>% slice(1:length(y_index_mis[[2]]))

fit_summary_mcelreath<-summary(fit.mcelreath.miss[[2]], probs=c(0.025,.5,.975))$summary

my_data <- as_tibble(fit_summary_mcelreath)
y_est_data_mcelreath<-my_data %>% slice(1:length(y_index_mis[[2]]))


y_combined<-cbind(y_est_data_mcelreath[,1],y_est_data_matt[,1])

##rename column names for clarity
colnames(y_combined)<-c("mcelreath", "matt")

##check that column names are correct

head(y_combined)

##Plot observed and estimated direct comparison
ggplot(data=y_combined, aes(x=matt, y=mcelreath))+
  geom_point( size=3)+
  geom_abline(intercept=0, slope=1)+
  theme_classic()+
  theme(legend.position="top")+
  ylab("Estimated GPP")+
  xlab("Observed GPP")+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))

