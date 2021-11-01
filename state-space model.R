# Load packages -----------------------------------------------------------

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


# Simulate data -----------------------------------------------------------

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

#Estimate 'real' light data
light<-lightest(time, 47.8762, -114.03, 105) #Flathead Bio Station just for fun

#Condense light to be relative to max light
light.rel<-light/max(light)


##GPP based on time-series model with known parameters
set.seed(550)
x<-NA
sdp <- 0.01
sdo<-0.1
phi <-0.8
b0<-0.1
b1<-0.1
x[1]<-.7

## Set the seed, so we can reproduce the results
set.seed(550)
## For-loop that simulates the state through time, using i instead of t,
for (t in 2:N){
  x[t] = b0+phi*x[t-1]+light.rel[t]*b1+rnorm(1, 0, sdp)
}

#Estimate observed data
set.seed(550)
sd.o <-rnorm(N, 0, sdo)
y <-x + sd.o

##Bind simulated data and time
y_full<-as.data.frame(cbind(x,y,time))

ggplot(data=y_full, aes(x=time, y=x))+
  geom_point(color="black", size=3)+
  geom_point(aes(y=y, x=time),color="red", size=3)+
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

# Model code for STAN -----------------------------------------------------
sink("ss.stan")

cat("
 

/*----------------------- Data --------------------------*/
  data {
    int N; // Length of state and observation time series
    vector[N] y; // Observations
    vector[N] light;
    }
  
/*----------------------- Parameters --------------------------*/
  
  parameters {
    vector[N] z; //latent state variable
    real<lower=0> sdp; // Standard deviation of the process equation
    real<lower=0> sdo; // Standard deviation of the observation equation
    real b0; // intercept
    real b1; // light coefficient 
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
      }
  

  /*----------------------- Model --------------------------*/
  model {
    // Prior distributions
    sdo ~ normal(0, 1);
    sdp ~ normal(0, 1);
    phi ~ beta(1,1);
    b0 ~ normal(0,5);
    b1 ~ normal(0,5);
    
    // Distribution for the first state
    z[1] ~ normal(y[1], 0.1);
  
    // Distributions for all other states
    for(t in 2:N){
       z[t] ~ normal(b0+ z[t-1]*phi+light[t]*b1, sdp);// process model with error
    }
       y ~ normal(exp(z), sdo); // observation model with fixed observation error
     }
  
  
    "
    ,fill=TRUE)
sink()
closeAllConnections()

#Prep model
model<-"ss.stan"
model<-stan_model(model)

#Create data object
data <- list(N = length(y),light=light.rel,y=y)

#Run Stan
fit<-rstan::sampling(object=model, data = data,  iter = 3000, chains = 4)#, control=list(max_treedepth = 15))


##Print and extract
fit_extract<-rstan::extract(fit)
print(fit, pars=c( "sdo","phi",  "sdp", "b0" ))

pairs(fit, pars=c("sdo","phi",  "sdp","b0"))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit)

##Traceplots

traceplot(fit, pars=c("phi"))

##Plot density plots
plot_sdo <- stan_dens(fit, pars="sdo") + geom_vline(xintercept =sdo)
plot_sdp <- stan_dens(fit, pars="sdp") + geom_vline(xintercept = sdp)
plot_phi <- stan_dens(fit, pars="phi") + geom_vline(xintercept = phi)
plot_b0 <- stan_dens(fit, pars="b0") + geom_vline(xintercept = b0)+xlab("Intercept")
plot_b1 <- stan_dens(fit, pars="b1") + geom_vline(xintercept = b1)+xlab("Light beta")
grid.arrange( plot_b0,plot_b1,plot_phi,plot_sdp,plot_sdo, nrow=2)
