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
set.seed(620)

###Simulate data
N<-366 #length of data
t<-0:N #time
time<-seq(as.POSIXct("2020/1/1"),as.POSIXct("2020/12/31"), by="day") #for light

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

##GPP based on time-series model with known paramters
x<-NA
x[1]=40 #for initialization
b0<-1 #intercept
phi<-0.8 #autoregressive parameter
b1<-1.5 #light parameter
sd.p<-1 #process error

#Estimate latent state
for (t in 2:N){
        x[t] = phi*x[t-1]+b1*light[t]+rnorm(1, 0, sd.p)
}
#Estimate observed data
sdo <-rnorm(N, 0, 1.5)
y <- x + sdo

sdo_sd<-mean(abs(sdo))

##Plot observed and latent
plot(1:N, x[1:366],
     pch=3, cex = 0.8, col="blue", ty="o", lty = 3,
     xlab = "t", ylab = expression(z[t]),
     xlim = c(0,N), ylim = c(min(y), max(y)+max(y)/5),
     las = 1)
points(1:N, y,
       pch = 19, cex = 0.7, col="red", ty = "o")
legend("top",
       legend = c("Obs.", "Latent states"),
       pch = c(3, 19),
       col = c("blue", "red"),
       lty = c(3, 1),
       horiz=TRUE, bty="n", cex=0.9)

##Force some missing data
#Note-the model will still work if no data is missing but the objects in the next block of code (i.e., "y_index_mis", etc.)
#still need to be created
y_miss<-y
#y_na<- which(y_miss %in% sample(y_miss,50)) ##the X in (y_miss,X) is the number of created missing data points 
#y_miss[y_na]<-NA

##Create objects with the location (index) of missing or observed data AND number
## of missing or observed data
y_index_mis <- which(is.na(y_miss)) # Identify the rows for each variable in airquality with NA
y_index_obs <- which(!is.na(y_miss)) # Identify the rows for each variable in airquality with observed data
y_nMiss <- length(y_index_mis)# How many NAs?
y_nObs <- length(y_index_obs) # How many NAs?

##Replace NAs with arbitrary number to make Stan happy
y_miss[is.na(y_miss)] <- -100

##Stan Code

sink("ss_reg_missing data.stan")

cat("
    data {
    int<lower=1> N; // Number of observations
    vector[N] y_miss; // Response variable, including missing values
    int y_nMiss; // number of missing values
    int y_index_mis[y_nMiss]; // index or location of missing values within the dataset
    real light[N]; //light data
    }
    
    parameters {
    vector[y_nMiss] y_imp;// Missing data
    vector[N] X; // latent state data
    real<lower = 0> sigma_proc; // process error
    real<lower = 0> sigma_obs; // observation error
    real<lower = 0, upper=1 > phi;  // auto-regressive parameter
    real b1; // light parameter 
    }
    
    transformed parameters { 
    vector[N] y;
    y=y_miss; // makes the data a transformed variable
    y[y_index_mis] =y_imp; // replaces missing data in y with estimated data
    } 
    
    model {
        X[1]~normal(y[1],sigma_proc); //set initial state
    
    for(i in 1:N){
    y[i] ~ normal(X[i], sigma_obs); // observation model
    }
    
    for (t in 2:N){
    X[t] ~ normal(phi*X[t-1]+b1*light[t], sigma_proc); // process model with unknown process error //regression model with AR errors
    
    }  
    
    // error priors
    sigma_proc ~ normal(0, 1); 
    sigma_obs ~ normal(0, 2);
    
    // single parameters priors 
    b1~normal(0, 5);
   
    //Prior for AR coefficient needs to be reparameterized
    phi~beta(1,1);
    }
 
generated quantities {
 vector[N] x_rep; // replications from posterior predictive dist
 vector [N] y_rep;
 x_rep[1]=normal_rng(y[1], 0.1);
 
 for (t in 2:N) {
 x_rep[t]=normal_rng(phi*X[t-1]+b1*light[t], sigma_proc);
 }
 for (t in 1:N) {
 y_rep[t]=normal_rng(x_rep[t], sigma_obs);
 }

  
 }
 
    "
,fill=TRUE)
sink()
closeAllConnections()

##Load data

data <- list(   N = N,
                y_nMiss = y_nMiss,
                y_nObs = y_nObs,
                y_index_mis =y_index_mis,
                y_index_obs = y_index_obs,
                y_miss= y_miss,
                light=light
                
)


##Run Stan
fit<- stan("ss_missing data.stan", data = data,  iter = 5000, chains = 4)

##Print and extract
fit_extract<-rstan::extract(fit)
print(fit, pars=c( "sigma_proc","phi", "b1", "sigma_obs" ))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit)

##Traceplots
traceplot(fit, pars=c( "sigma_proc", "sigma_obs" ))
traceplot(fit, pars=c("phi", "b1"))

##Posterior densities comopared to known parameters
plot_sdo <- stan_dens(fit, pars="sigma_obs") + geom_vline(xintercept =1.5)
plot_sdp <- stan_dens(fit, pars="sigma_proc") + geom_vline(xintercept = sd.p)
plot_phi <- stan_dens(fit, pars="phi") + geom_vline(xintercept = phi)
plot_b1 <- stan_dens(fit, pars="b1") + geom_vline(xintercept = b1)
grid.arrange(plot_b1, plot_phi,plot_sdp,plot_sdo, nrow=2)

##Posterioir Predictive check and test statistic
yrep1 <- fit_extract$y_rep
samp100 <- sample(nrow(yrep1), 100)
ppc_dens_overlay(y_miss, yrep1[samp100, ]) 
ppc_stat(y_miss, yrep1[samp100, ])

##Shineystan for everything else
launch_shinystan(fit)



######Compare 'missing' data with simulated data
##Create object with estimated missing data

fit_summary<-summary(fit, probs=c(0.05,.5,.95))$summary

my_data <- as_tibble(fit_summary)
y_est_data<-my_data %>% slice(1:length(y_index_mis))

##Create object with observed data

y_miss_index<-data$y_index_mis
date<-1:TT
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
  ylab("Estimated Temperature")+
  xlab("Observed Temperature")+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))

##Plot observed and estimated time series
time2=2:TT
raw.data=as.data.frame(cbind(time2, y))

ggplot(data=y_combined, aes(x=date, y=estimated))+
  geom_line(data=raw.data, aes(x=time2, y=y, color="black"), size=1)+
  geom_point(aes(color="red"), size=3)+
  geom_errorbar(data=y_combined,aes(x=date,ymin=low, ymax=high))+ 
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  scale_color_identity(guide = "legend", name=element_blank(), labels=c("Raw data", "Estimated points"))+
  theme(legend.position="top")+
  ylab("Sim data")





