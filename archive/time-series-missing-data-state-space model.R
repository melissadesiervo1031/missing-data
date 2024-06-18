library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(dplyr)
library(shinystan)
library(faux)
library(bayesplot)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(data.table)

###Simulate data
N<-731 #length of data
t<-0:N #time
time<-seq(as.POSIXct("2020/1/1"),as.POSIXct("2021/12/31"), by="day") #for light

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
light<-log(light)

##Turbidity--Duwamish River, Tukwila, Washington (https://www.sciencebase.gov/catalog/item/5a1dba7de4b09fc93dd7c022)
fake.turb<-read.csv("fake_turb.csv")
length(fake.turb$turb)
turb<-fake.turb$turb

##GPP based on time-series model with known paramters
set.seed(360)
x<-NA
b0<-1#intercept
phi<-0.9 #autoregressive parameter
b1<-0.8 #light parameter
#b2<--.05 #turbidity parameter
sd.p<-0.5#process error
x[1]<-25 #for initialization

#Estimate latent state
for (t in 2:N){
       x[t] = b0+phi*x[t-1]+b1*light[t]+rnorm(1, 0, sd.p)
}
#Estimate latent state
#for (t in 2:N){
#        x[t] = b0+phi*x[t-1]+rnorm(1, 0, sd.p)
#}


#Estimate observed data
set.seed(360)
sd.o<-0.5
sdo <-rnorm(N, 0, sd.o)
y <-x + sdo

sdo_sd<-mean(abs(sdo))

##Plot observed and latent
plot(1:N, x[1:N],
     pch=3, cex = 0.8, col="blue", ty="o", lty = 3,
     xlab = "t", ylab = expression(z[t]),
     xlim = c(0,N), ylim = c(min(y), max(y)),
     las = 1)
points(1:N, y,
       pch = 19, cex = 0.7, col="red", ty = "o")
legend("top",
       legend = c("Obs.", "Latent states"),
       pch = c(3, 19),
       col = c("blue", "red"),
       lty = c(3, 1),
       horiz=TRUE, bty="n", cex=0.9)

#Center data
y<-(y-mean(y))/sd(y)
mean(y)
sd(y)

x<-(x-mean(x))/sd(x)
mean(x)
sd(x)


plot(1:N,y)
lines(1:N,x, color="blue")
##Force some missing data
#Note-the model will still work if no data is missing but the objects in the next block of code (i.e., "y_index_mis", etc.)
#still need to be created

set.seed(430)
#Store simulated data in a list 
y_miss<-list(y,y,y,y,y,y,y,y)

#vector of missing number amounts
missing_n<-c(0,10,25,50,75,100,150,175)

#Create a list to store missing data integers
y_missing_integers<-c()

#Create initial missing data integer set. Then build on with the for loop
y_missing_integers[[1]]<- which(y %in% sample(y,missing_n[1]))

#Adds the previous missing data integers to the next set. This makes them nested.
#i.e, all the missing data 10 integers are added to 15 more to make missing data 15.
for(i in 2:length(missing_n)){
z<-which(y %in% sample(y,(missing_n[i]-missing_n[i-1]), replace = FALSE))
y_missing_integers[[i]]<-c(x,y_missing_integers[[i-1]])
}

#Assign NAs to the missing data in each list
for(i in 1:length(missing_n)){
  z<-y_miss[[i]]
  z[y_missing_integers[[i]]]<-NA
  z[1]<-y[1]
  y_miss[[i]]<-z
}


#Check that each list has the right amount of missing data (or close to it...)
miss<-sapply(y_miss, function(y) sum(length(which(is.na(y)))))
prop.miss<-round(miss/length(y)*100)
prop.miss

##Create objects with the location (index) of missing or observed data AND number
## of missing or observed data
#FOR SINGLE VECTOR
y_index_mis <- which(is.na(y_miss)) # Identify the rows for each variable in airquality with NA
y_index_obs <- which(!is.na(y_miss)) # Identify the rows for each variable in airquality with observed data
y_nMiss <- length(y_index_mis)# How many NAs?
y_nObs <- length(y_index_obs) # How many NAs?

##Create objects with the location (index) of missing or observed data AND number
## of missing or observed data
#FOR LIST OF VECTORS
y_index_mis <- lapply(y_miss,function(var){which(is.na(var))}) # Identify the rows for each variable in airquality with NA
y_index_obs <- lapply(y_miss,function(var){which(!is.na(var))}) # Identify the rows for each variable in airquality with observed data
y_nMiss <- lapply(y_index_mis,function(var){length(var)}) # How many NAs?
y_nObs <- lapply(y_index_obs,function(var){length(var)}) # How many NAs?


##Replace NAs with arbitrary number to make Stan happy
for(i in 1:length(missing_n)){
  x<-y_miss[[i]]
  x[y_missing_integers[[i]]]<--100 #arbitrary number
  y_miss[[i]]<-x
} 

##Stan Code

sink("ss_missing data.stan")

cat("
    data {
    int<lower=1> N; // Number of observations
    vector[N] y_miss; // Response variable, including missing values
    int y_nmiss; // number of missing values
    int y_ind_mis[y_nmiss]; // index or location of missing values within the dataset
    real light[N]; //light data
    //real turb[N]; //turb data
    }
    
    parameters {
    vector[y_nmiss] y_imp;// Missing data
    vector[N] X; // latent state data
    real<lower = 0> sigma_proc; // process error
    real<lower = 0> sigma_obs; // observation error
    real<lower = 0, upper=1 > phi;  // auto-regressive parameter
    real b0; // intercept
    real b1; // light parameter
    //real b2; // turb paramater
    }
    
    transformed parameters { 
    vector[N] y;
    y=y_miss; // makes the data a transformed variable
    y[y_ind_mis] =y_imp; // replaces missing data in y with estimated data
    } 
    
    model {
    X[1]~normal(y[1],sigma_proc); //set initial state
    
   
    y ~ normal(X, sigma_obs); // observation model
    
    
    for (t in 2:N){
    X[t] ~ normal(b0+phi*X[t-1]+b1*light[t], sigma_proc); // process model with unknown process error //regression model with AR errors
    
    }  
    
    // error priors
    sigma_proc ~ normal(0, 1); 
    sigma_obs ~ normal(0, 1);
    
    // single parameters priors 
    b0~normal(0, 5);
    b1~normal(0, 5);
  
    //Prior for AR coefficient needs to be reparameterized
    phi~beta(1,1);
    }
 
generated quantities {
 vector[N] x_rep; // replications from posterior predictive dist
 vector [N] y_rep;
 x_rep[1]=normal_rng(y[1], 0.1);
 
 for (t in 2:N) {
 x_rep[t]=normal_rng(b0+phi*X[t-1]+b1*light[t], sigma_proc);
 }
 for (t in 1:N) {
 y_rep[t]=normal_rng(x_rep[t], sigma_obs);
 }

  
 }
 
    "
,fill=TRUE)
sink()
closeAllConnections()

fit.miss <- vector("list",length(missing_n))

for(i in 1:1){
##Load data
data <- list(   N = N,
                y_nmiss = unlist(y_nMiss[i]),
                y_nobs = unlist(y_nObs[i]),
                y_ind_mis = unlist(y_index_mis[i]),
                y_ind_obs = unlist(y_index_obs[i]),
                y= unlist(y_miss[i]),
                light=light)

##Run Stan
fit.miss[[1]]<- stan("ss.stan", data = data,  iter = 4000, chains = 4, control=list(max_treedepth = 15))
}

##Pull param estimates into list
fit_summary_pars <- vector("list",length(missing_n))
for (i in 1:length(missing_n)){
   fit_summary_pars[[i]]<-(summary(fit.miss[[i]], pars=c("sigma_proc","phi", "b1", "sigma_obs"), probs=c(0.05,.5,.95))$summary)
   }
fit_summary_pars[[1]]

##Unlist,cleanup, and add factors
fit_summary<-as.data.frame(do.call(rbind.data.frame, fit_summary_pars)) #Unlist
fit_summary$param<-rep(c("sigma_proc","phi", "b1", "sigma_obs"), times=length(missing_n) )#add parameter factor
fit_summary$prop.missing<-rep(prop.miss, times=c(4,4,4,4,4,4,4,4)) #add prop of missing data
row.names(fit_summary)<-1:(length(missing_n)*4) #remove row names
fit_summary$param<-as.factor(fit_summary$param) #make factor
fit_summary$prop.missing<-as.factor(fit_summary$prop.missing) #make factor
summary(fit_summary) #check it looks good

known.data<-c(b1,phi,sd.o,sd.p)
known<-as.data.frame(known.data)

known.param<-c("b1","phi",  "sigma_obs", "sigma_proc")
#known.missing<-c(5, 5, 5,5)
known<-as.data.frame(cbind(known.data, known.param))
known$known.data<-as.numeric(as.character(known$known.data))

ggplot(data=fit_summary, aes(x=mean, y=prop.missing ))+
  geom_point(aes( color=param, group=param),size=3,position=position_dodge(0.5))+
  theme_classic()+
  geom_errorbar(aes(xmin=mean-sd, xmax=mean+sd,color=param, group=param), width=0.2, size=0.5,position=position_dodge(0.5))+
  theme(legend.position="top")+
  ylab("Proportion of Missing Data")+
  xlab("Parameter Estimate")+
  scale_color_manual(values=c("blue", "green", "orange", "black"))+
  scale_x_continuous(limits=c(0.2,1))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  geom_vline(xintercept = known.data,color=c("blue", "green", "orange", "black"))
  


write.csv(fit_summary_pars_10, "Rstan_missing data 10_proc_0.5.csv")
##Print and extract
fit_extract<-rstan::extract(fit[[1]])
print(fit[[2]], pars=c( "sigma_proc","phi", "b1", "sigma_obs" ))

pairs(fit, pars=c("sigma_proc","phi", "b1", "sigma_obs"))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit.miss[[1]])

##Traceplots

traceplot(fit.miss[[1]], pars=c("phi", "b1","sigma_proc", "sigma_obs"))

##Posterior densities comopared to known parameters
plot_sdo <- stan_dens(fit.miss[[1]], pars="sigma_obs") + geom_vline(xintercept =sd.o)
plot_sdp <- stan_dens(fit.miss[[1]], pars="sigma_proc") + geom_vline(xintercept = sd.p)
plot_phi <- stan_dens(fit.miss[[1]], pars="phi") + geom_vline(xintercept = phi)
plot_b1 <- stan_dens(fit.miss[[1]], pars="b1") + geom_vline(xintercept = b1)+xlab("Light beta")
grid.arrange(plot_b1, plot_phi,plot_sdp,plot_sdo, nrow=2)

##Posterioir Predictive check and test statistic
yrep1 <- fit_extract$y_rep
samp100 <- sample(nrow(yrep1), 100)
ppc_dens_overlay(y, yrep1[samp100, ]) 
ppc_stat(y, yrep1[samp100, ])

##Shineystan for everything else
launch_shinystan(fit)

######Compare 'missing' data with simulated data
##Create object with estimated missing data

fit_summary<-summary(fit.10, probs=c(0.05,.5,.95))$summary

my_data <- as_tibble(fit_summary)
y_est_data<-my_data %>% slice(1:length(y_index_mis))

##Create object with observed data

y_miss_index<-data$y_index_mis
date<-1:N
y_obs_data<-as.data.frame(cbind(date[data$y_index_mis],y[data$y_index_mis]))
colnames(y_obs_data)<-c("date","y")
#y_obs_data$date<-as.Date(y_obs_data$date, origin = "1969-12-30")
#x_data<-summary(fit_extract$x, probs=c(0.05,.5,.95))
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
  ylab("Estimated GPP")+
  xlab("Observed GPP")+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))


###############################
#####Multiple Imputations######
###############################
ts<-as.Date(time,format="%m/%d/%Y")

data.amelia<- vector("list",length(missing_n))
for(i in 1:length(missing_n)){
y_miss1<-na_if(y_miss[[i]], -100)
y_miss1[1]<-y[1]
z<-as.data.frame(cbind(y_miss1, ts))
z2<-amelia(z, m = 5, p2s=1, ts="ts", polytime = 2, emperi=.1*nrow(ts))
summary(z2)
data.amelia[[i]]<-z2
}
str(z2$imputations,1)

fit.amelia <- vector("list",length(missing_n))

for(i in 2:length(missing_n)){##two because the list holding the lists starts with a list of misc data that arent the estimates
  for (g in 1:5){
  ##Load data
  data <- list(   N = N,
                  y_miss= data.amelia[[2]]$imputations[[1]]$y_miss1,
                  light=light)

}
}
  ##Run Stan
  fit.1[[i]]<- stan("ss_missing data.stan", data = data,  iter = 4000, chains = 4, control=list(max_treedepth = 15))


##Stan Code

sink("ss.stan")

cat("
    data {
    int<lower=1> N; // Number of observations
    vector[N] y; // Response variable, including missing values
    vector[N] light; //light data
   
    }
    
    parameters {
    vector[N] X; // latent state data
    real<lower = 0> sigma_proc; // process error
    real<lower = 0> sigma_obs; // observation error
    real<lower = 0, upper=1 > phi;  // auto-regressive parameter
    real b0; // intercept
    real b1; // light parameter
      }
    
    model {
    X[1]~normal(y[1],0.01); //set initial state
    
   
    y ~ normal(X, sigma_obs); // observation model
    
    
    for (t in 2:N){
    X[t] ~ normal(b0+phi*X[t-1]+b1*light[t], sigma_proc); // process model with unknown process error 
    
    }  
    
    // error priors
    sigma_proc ~ normal(0, 1); 
    sigma_obs ~ normal(0, 1);
    
    // single parameters priors 
    b0~normal(0, 5);
    b1~normal(0, 5);
  
    //Prior for AR coefficient needs to be reparameterized
    phi~beta(1,1);
    }
 
generated quantities {
 vector[N] x_rep; // replications from posterior predictive dist
 vector [N] y_rep;
 x_rep[1]=normal_rng(y[1], 0.1);
 
  x_rep[1]=normal_rng(y[1], 0.1);
 y_rep[1]=normal_rng(y[1], 0.1);
 
 for (t in 2:N) {
 x_rep[t]=normal_rng(b0+phi*X[t-1]+b1*light[t], sigma_proc);
 y_rep[t]=normal_rng(x_rep[t], sigma_obs);
 }
 }
 
    "
    ,fill=TRUE)
sink()
closeAllConnections()

#for(i in 1:1){
  ##Load data
  data <- list(   N = N,
                  y= y,
                  light=light)
  
  ##Run Stan
  fit.miss<- stan("ss.stan", data = data,  iter = 1000, chains = 4)
  ##Traceplots
  
  traceplot(fit.miss, pars=c("b0","b1", "phi", "sigma_proc", "sigma_obs"))
  
  ##Posterior densities comopared to known parameters
  plot_sdo <- stan_dens(fit.miss, pars="sigma_obs") + geom_vline(xintercept =sd.o)
  plot_sdp <- stan_dens(fit.miss, pars="sigma_proc") + geom_vline(xintercept = sd.p)
  plot_phi <- stan_dens(fit.miss, pars="phi") + geom_vline(xintercept = phi)
  plot_b1 <- stan_dens(fit.miss, pars="b1") + geom_vline(xintercept = b1)+xlab("Light beta")
  plot_b0 <- stan_dens(fit.miss, pars="b0") + geom_vline(xintercept = b0)+xlab("Intercept")
  grid.arrange( plot_b0, plot_phi,plot_sdp,plot_sdo, nrow=2)
  
