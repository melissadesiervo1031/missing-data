#Load packages

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
library(R.utils)
library(dlm)

# PHI + SDP (process model)---------------------------------------------------------------------

#Simulate data

N<-365 #length of data
t<-1:N #time
time<-seq(as.POSIXct("2020/1/1"),as.POSIXct("2020/12/30"), by="day") #for light

##GPP based on time-series model with known parameters
set.seed(550)
z<-NA
sdp <- 0.01
sdo<-0.01
phi <-0.8
z[1]<-0

## Set the seed, so we can reproduce the results
set.seed(550)
## For-loop that simulates the state through time, using i instead of t,
for (t in 2:N){
  z[t] = phi*z[t-1]+rnorm(1, 0, sdp)
}

y_full<-data.frame(z,t)

ggplot(data=y_full, aes(x=time, y=z))+
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

#STAN model
sink("ss1.stan")

cat("
 

/*----------------------- Data --------------------------*/
  data {
    int N; // Length of state and observation time series
    vector[N] z; // Observations
    real z0;
    }
  
/*----------------------- Parameters --------------------------*/
  
  parameters {
    real<lower=0> sdp; // Standard deviation of the process equation
    real b0;
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
      }
  

  /*----------------------- Model --------------------------*/
  model {
    // Prior distributions
    sdp ~ normal(0, 1);
    phi ~ beta(1,1);
    b0 ~ normal(0,5);
    
    // Distribution for the first state
    z[1] ~ normal(z0, sdp);
  
    // Distributions for all other states
    for(t in 2:N){
       z[t] ~ normal(b0+phi*z[t-1], sdp);// process model with error
    }
  }
    
    "
  
    ,fill=TRUE)
sink()
closeAllConnections()

#Prep model
model<-"ss1.stan"
model<-stan_model(model)

#Create data object
data <- list(N = length(z),z=z, z0=z[1])

#Run Stan
fit<-rstan::sampling(object=model, data = data,  iter = 3000, chains = 4)#, control=list(max_treedepth = 15))


saveRDS(fit, file = "process.RDS") 

##Print and extract
fit_extract<-rstan::extract(fit)
print(fit, pars=c( "phi","sdp"))

pairs(fit, pars=c("phi","sdp"))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit)

##Traceplots

traceplot(fit, pars=c("phi", "sdp", "b0"))

##Density plots
plot_sdp <- stan_dens(fit, pars="sdp") + geom_vline(xintercept = sdp)
plot_phi <- stan_dens(fit, pars="phi") + geom_vline(xintercept = phi)
plot_b0 <-stan_dens(fit, pars="b0") + geom_vline(xintercept = 0)
grid.arrange(plot_phi,plot_sdp, plot_b0, nrow=2)


#SDP + SDO AM apporach  (state-space model)-------------------------------------------------


#Simulate data

N<-1000#length of data
z<-numeric(N+1) 
t<-numeric(N+1) 
t<-0:500
##GPP based on time-series model with known parameters
set.seed(550)
sdp <- 0.1
sdo<-0.01

phi<-0.8

## Set the seed, so we can reproduce the results
set.seed(550)
sd.p<-rnorm(N, 0, sdp)
mean(sd.p)
sd(sd.p)
## For-loop that simulates the state through time
for (i in 1:N){
  z[i+1] =z[i]*phi+sd.p[i]
}

#Estimate observed data
set.seed(553)
sd.o <-rnorm(N, 0, sdo)
mean(sd.o)
sd(sd.o)
y <-z[2:(N+1)]+ sd.o



#y<-insert(y, 1, NA)
#sd.o<-insert(sd.o, 1, NA)
#sd.p<-insert(sd.p, 1, NA)
##Bind simulated data and time
y_full_am<-as.data.frame(cbind(t,z,y, sd.o, sd.p))


##Plot
plot(1:length(y), y,
     pch=3, cex = 0.8, col="blue", ty="o", lty = 3,
     xlab = "t", ylab = expression(y[t]),
     xlim = c(0,N), #ylim = c(min(y), max(y+max(y)/5)),
     las = 1)
points(0:N, z,
       pch = 19, cex = 0.7, col="red", ty = "o")
legend("top",
       legend = c("Obs.", "True states"),
       pch = c(3, 19),
       col = c("blue", "red"),
       lty = c(3, 1),
       horiz=TRUE, bty="n", cex=0.9)

#Stan model
sink("ss_2.stan")

cat("
 

/*----------------------- Data --------------------------*/
  data {
    int N; // Length of state and observation time series
    vector[N] y; // Observations
    real z0; //intital mean of prior state value
    

    }
  
/*----------------------- Parameters --------------------------*/
  
  parameters {
    
    vector[N] z; //latent state variable
    real<lower =0> sdp; // Standard deviation of the process equation
    real<lower =0> sdo; // Standard deviation of the process equation
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    real b0; 
    
    }
  

  /*----------------------- Model --------------------------*/
  model {
    // Prior distributions
    phi ~ beta(1,1);
    b0 ~ normal(0, 5);
    sdp ~normal(0,1);
    sdo ~normal(0.01,0.001);
    
    // Observation error model       
 
   for(i in 1:N){
   y[i] ~ normal(z[i], sdo); // observation model with fixed observation error
   }
  
    
    // Distribution for the zero and first state
    
    z[1] ~ normal(z0,sdp);
   
    // Distributions for states
    for(i in 2:N){
       z[i] ~ normal(z[i-1]*phi+b0, sdp);// process model with error
        }
   
 
 }
"
    ,fill=TRUE)
sink()
closeAllConnections()


#Prep model
model<-"ss_2.stan"
model<-stan_model(model)

#Create data object
data <- list(N = length(y),y=y, z0=y[1])

#Run Stan
fit.2<-rstan::sampling(object=model, data = data,  iter = 3000, chains = 4)#, control=list(max_treedepth = 15))

fit_extract<-rstan::extract(fit)
print(fit.2, pars=c( "sdo", "sdp", "phi","b0"))

##Check pairs plot
pairs(fit, pars=c("sdo", "sdp", "phi", "b0"))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit)

##Traceplots

traceplot(fit.2, pars=c("sdo","sdp", "phi", "b0"))

##Plot density plots
plot_sdo <- stan_dens(fit.2, pars="sdo") + geom_vline(xintercept =sdo)
plot_sdp <- stan_dens(fit.2, pars="sdp") + geom_vline(xintercept = sdp)
plot_phi <- stan_dens(fit.2, pars="phi") + geom_vline(xintercept = phi)
plot_b0 <- stan_dens(fit.2, pars="b0") + geom_vline(xintercept =0)
grid.arrange(plot_sdp, plot_sdo, plot_phi, plot_b0, nrow=2)

saveRDS(fit.2, file = "ss_AM_fixedsdo.RDS") 

#SDP + SDO-Hawari (state-space model)-------------------------------------------------


#Simulate data

N<-500#length of data
z<-numeric(N+1) 
t<-numeric(N+1) 
t<-0:500
##GPP based on time-series model with known parameters
set.seed(553)
sdp <- 0.01
sdo<-0.001
phi<-.8

## Set the seed, so we can reproduce the results
set.seed(550)
sd.p<-rnorm(N, 0, sdp)
## For-loop that simulates the state through time
for (i in 1:N){
  z[i+1] =z[i]*phi+sd.p[i]
}

#Estimate observed data
set.seed(500)
sd.o <-rnorm(N, 0, sdo)
y <-z[2:(N+1)]+ sd.o



#y<-insert(y, 1, NA)
#sd.o<-insert(sd.o, 1, NA)
#sd.p<-insert(sd.p, 1, NA)
##Bind simulated data and time
y_full_am<-as.data.frame(cbind(t,z,y, sd.o, sd.p))


##Plot
plot(1:length(y), y,
     pch=3, cex = 0.8, col="blue", ty="o", lty = 3,
     xlab = "t", ylab = expression(y1[t]),
     xlim = c(0,N), #ylim = c(min(y), max(y+max(y)/5)),
     las = 1)
points(0:N, z,
       pch = 19, cex = 0.7, col="red", ty = "o")
legend("top",
       legend = c("Obs.", "True states"),
       pch = c(3, 19),
       col = c("blue", "red"),
       lty = c(3, 1),
       horiz=TRUE, bty="n", cex=0.9)

#Stan model
sink("ss_3.stan")

cat("
 

/*----------------------- Data --------------------------*/
  data {
    int<lower=1> N; // Length of state and observation time series
    vector[N] y; // Observations
    real m0; //intital mean of prior state value
    cov_matrix[1] sdp0; //variance of prior distribution
    real n0; //intital mean of prior state value
    cov_matrix[1] sdo0; //variance of prior distribution
    }
  
/*----------------------- Parameters --------------------------*/
  
  parameters {
    real z0; //state [0];
    real y0; //state [0];
    vector[N] z; //latent state variable
    cov_matrix[1] sdp; // Standard deviation of the process equation
    cov_matrix[1] sdo; // Standard deviation of the process equation
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    real b0; 
    
    }
  

  /*----------------------- Model --------------------------*/
  model {
    // Prior distributions
    phi ~ beta(1,1);
    b0 ~ normal(0, 5);
   
    
     // Observation error model       
    y0 ~ normal(n0,sdo0[1,1]);
    y[1] ~ normal(y0,sdo[1,1]);
    
    for(i in 2:N){
    y[i] ~ normal(z[i], sdo[1,1]); // observation model with fixed observation error
    }
    
    // Distribution for the zero and first state
    z0 ~ normal(m0,sdp0[1,1]);
    z[1] ~ normal(z0,sdp[1,1]);
   
    // Distributions for states
    for(i in 2:N){
       z[i] ~ normal(z[i-1]*phi+b0, sdp[1,1]);// process model with error
        }
   
 
 }
"
    ,fill=TRUE)
sink()
closeAllConnections()

#Prep model
model<-"ss_3.stan"
model<-stan_model(model)

#Create data object
data <- list(N=N, y=y, m0=0.01, sdp0=matrix(0.1, nrow=1, ncol=1),n0=0, sdo0=matrix(1, nrow=1, ncol=1))


#Run Stan
fit<-rstan::sampling(object=model, data = data,  iter = 3000, chains = 4)#, control=list(max_treedepth = 15))

#saveRDS(fit, file = "ss_AM.RDS") 

##Print and extract
fit_extract<-rstan::extract(fit)
print(fit, pars=c( "sdo", "sdp", "phi","b0"))

##Check pairs plot
pairs(fit, pars=c("sdo", "sdp", "phi", "b0"))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit)

##Traceplots

traceplot(fit, pars=c("sdo","sdp", "phi", "b0"))

##Plot density plots
plot_sdo <- stan_dens(fit, pars="sdo") + geom_vline(xintercept =sdo)
plot_sdp <- stan_dens(fit, pars="sdp") + geom_vline(xintercept = sdp)
plot_phi <- stan_dens(fit, pars="phi") + geom_vline(xintercept = phi)
plot_b0 <- stan_dens(fit, pars="b0") + geom_vline(xintercept =0)
grid.arrange(plot_sdp,plot_sdo, plot_phi, plot_b0, nrow=2)




##Create object with estimated missing data

fit_summary<-summary(fit, probs=c(0.05,.5,.95))$summary
my_data <- as_tibble(fit_summary)
y_est_data<-my_data[2:501,]

plot(y_est_data$mean~x)

##Merge observed and estimated

y_combined<-cbind(y_full_am$t[2:501],y_full_am$z[2:501],y_full_am$y[1:500],y_est_data)

##rename column names for clarity

colnames(y_combined)<-c("time","obs.z", "obs.y","est.z", "sd.z", "low", "median","high", "n_eff", "rhat")

##check that column names are correct

head(y_combined)

##Plot observed and estimated direct comparison

ggplot(data=y_combined, aes(x=obs.z, y=est.z))+
  geom_point( size=3)+
  geom_abline(intercept=0, slope=1)+
  theme_classic()+
  #geom_errorbar(aes(ymin=low, ymax=high), width=0.5)+
  theme(legend.position="top")+
  ylab("Estimated state")+
  xlab("Observed state")+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))




# SDP + SDO MT/ROH approach  (state-space model)-------------------------------------------------


#Simulate data

N<-500 #length of data
t<-1:500 #time


##GPP based on time-series model with known parameters
set.seed(550)
z<-NA
sdp <- 0.01
sdo<-0.01
z[1]<-0

## Set the seed, so we can reproduce the results
set.seed(550)
sd.p<-rnorm(N, 0, sdp)
## For-loop that simulates the state through time
for (i in 2:N){
  z[i] =z[i-1]+sd.p[i]
}

#Estimate observed data
set.seed(550)
sd.o <-rnorm(N, 0, sdo)
y <-z+ sd.o

##Bind simulated data and time
y_full<-as.data.frame(cbind(t,z,y,sd.o,sd.p))

ggplot(data=y_full, aes(x=t, y=z))+
  geom_point(color="black", size=3)+
  geom_point(aes(y=y, x=t),color="red", size=3)+
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

plot(z,y)

#Stan model
sink("ss_1.stan")

cat("
 

/*----------------------- Data --------------------------*/
  data {
    int N; // Length of state and observation time series
    vector[N] y; // Observations
    real z0; // Mean of prior distribution
    cov_matrix[1] cz0; // Variance of prior distribution
    }
  
/*----------------------- Parameters --------------------------*/
  
  parameters {
    vector[N] z; //latent state variable
    real<lower=0> sdp; // Standard deviation of the process equation
    real<lower=0> sdo; // Standard deviation of the observation equation
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    real b0; 
      }
  

  /*----------------------- Model --------------------------*/
  model {
    // Prior distributions
    sdo ~ normal(0, 1);
    sdp ~ normal(0, 1);
    phi ~ beta(1,1);
    b0  ~ normal(0,5);
    
    // Distribution for the first state
    z[1] ~ normal(z0,sdp);
  
    // Distributions for all other states
    for(i in 2:N){
       z[i] ~ normal(z[i-1]*phi+b0, sdp);// process model with error
    }
    
    for(i in 1:N){
       y[i] ~ normal(z[i], sdo); // observation model with fixed observation error
    }

 }
"
    ,fill=TRUE)
sink()
closeAllConnections()

#Prep model
model<-"ss_1.stan"
model<-stan_model(model)

#Create data object
data <- list(N = length(y),y=y, z0=z[1])

#Run Stan
fit.2<-rstan::sampling(object=model, data = data,  iter = 3000, chains = 4)#, control=list(max_treedepth = 15))

saveRDS(fit.2, file = "state_space_3.RDS") 
##Print and extract
fit_extract<-rstan::extract(fit)
print(fit, pars=c( "sdo","sdp", "phi", "b0"))

pairs(fit, pars=c("sdo","sdp"))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit)

##Traceplots

traceplot(fit.2, pars=c("sdo", "sdp"))

##Plot density plots
plot_sdo <- stan_dens(fit, pars="sdo") + geom_vline(xintercept =sdo)
plot_sdp <- stan_dens(fit, pars="sdp") + geom_vline(xintercept = sdp)
plot_phi <- stan_dens(fit, pars="phi") + geom_vline(xintercept = phi)
plot_b0 <- stan_dens(fit, pars="b0") + geom_vline(xintercept =0)
grid.arrange(plot_sdp,plot_sdo,plot_phi, plot_b0,nrow=2)


# FULL MODEL --------------------------------------------------------------



#Simulate data

N<-365 #length of data
z<-numeric(N+1) 
tt<-0:N
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
sdp <- 0.1
sdo<-0.01
phi <-0.8
b1<-0.6
b0<-2

z[1]<-10
## Set the seed, so we can reproduce the results
set.seed(5)
## For-loop that simulates the state through time, using i instead of t,
for (i in 1:365){
  z[i+1] = b0+z[i]*phi+light.rel[i]*b1+rnorm(1, 0, sdp)
}

#Estimate observed data
set.seed(4)
sd.o <-rnorm(N, 0, sdo)
y <-z[2:(N+1)] + sd.o
y[N+1]<-NA

##Bind simulated data and time
y_full<-as.data.frame(cbind(z,y,tt))
y_full$tt<-as.factor(y_full$tt)

ggplot(data=y_full, aes(x=tt, y=z))+
  geom_point(aes(color="state"), size=2)+
  geom_point(aes(y=y, x=tt,color="observed"), size=2)+
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
  scale_color_manual(name=c("observed","state"), values = c("red","black"))
  



#Stan model
sink("ss_light.stan")

cat("
 

/*----------------------- Data --------------------------*/
  data {
    int N; // Length of state and observation time series
    vector[N] y; // Observations
    vector[N] light;
    real z0;
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
    phi ~ beta(1,1);
    b0 ~ normal(0, 5);
    b1 ~ normal(0,5);
    sdp ~normal(0,1);
    sdo ~normal(0.01,0.001);
    
    // Observation error model       
 
   for(i in 1:N){
   y[i] ~ normal(z[i], sdo); // observation model with fixed observation error
   }
  
    
    // Distribution for the zero and first state
    
    z[1] ~ normal(z0,sdp);
   
    // Distributions for states
    for(i in 2:N){
       z[i] ~ normal(b0+ z[i-1]*phi+light[i]*b1, sdp);// process model with error
        }
   
 
 }
  
      "
    
    ,fill=TRUE)
sink()
closeAllConnections()

#Prep model
model<-"ss_light.stan"
model<-stan_model(model)

#Create data object
data <- list(N = length(y[1:365]),y=y[1:365], z0=y[1], light=light.rel)

#Run Stan
fit<-rstan::sampling(object=model, data = data,  iter = 3000, chains = 4)#, control=list(max_treedepth = 15))


##Print and extract
fit_extract<-rstan::extract(fit)
print(fit, pars=c( "sdo","phi", "sdp", "b0","b1" ))

pairs(fit, pars=c("sdo","phi", "sdp","b0","b1","lp__"))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit)

##Traceplots

traceplot(fit, pars=c("sdo","phi", "sdp", "b1", "b0"))

##Plot density plots
plot_sdo <- stan_dens(fit, pars="sdo") + geom_vline(xintercept =sdo)
plot_sdp <- stan_dens(fit, pars="sdp") + geom_vline(xintercept = sdp)
plot_phi <- stan_dens(fit, pars="phi") + geom_vline(xintercept = phi)
plot_b0 <- stan_dens(fit, pars="b0") + geom_vline(xintercept = 2)+xlab("Intercept")
plot_b1 <- stan_dens(fit, pars="b1") + geom_vline(xintercept = b1)+xlab("Light beta")
grid.arrange(plot_b0,plot_b1,plot_phi,plot_sdp,plot_sdo, nrow=2)

##Create object with estimated missing data

fit_summary<-summary(fit, probs=c(0.05,.5,.95))$summary
my_data <- as_tibble(fit_summary)
y_est_data<-my_data[1:365,]


##Merge observed and estimated

y_combined<-cbind(y_full$z[2:366],y_full$y[1:365],y_full$tt[1:365],y_est_data)

##rename column names for clarity

colnames(y_combined)<-c("obs.z", "obs.y", "time", "est.z", "se.z","sd.z", "low", "median","high", "n_eff", "rhat")

##check that column names are correct

head(y_combined)

##Plot observed and estimated direct comparison

ggplot(data=y_combined, aes(x=obs.z, y=est.z))+
  geom_point( size=3)+
  geom_abline(intercept=0, slope=1)+
  theme_classic()+
  #geom_errorbar(aes(ymin=low, ymax=high), width=0.5)+
  theme(legend.position="top")+
  ylab("Estimated state")+
  xlab("Observed state")+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))






