# setwd('C:/Users/alice.carter/git/missing-data/')

#Load packages

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(bayesplot)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(data.table)

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

# center at expected value:
EV = (b0*(1-phi) + b1*light.rel[1])/(1-phi)
z[1]<-EV
## Set the seed, so we can reproduce the results
# set.seed(5)
## For-loop that simulates the state through time, using i instead of t,
# In this process equation, I have centered the data to help remove the equifinality
# between b0 and phi. See:
# https://stats.stackexchange.com/questions/348430/what-is-the-reason-for-not-including-an-intercept-term-in-ar-and-arma-models
# with no intercept: y_t = phi * y_t-1 + b1*light_t + w_t
# center y by subtracting the expected intercept (b0):
#                   y_t - b0 = phi*(y_t-1 - b0) + b1 * light_t + w_t
# rearrange:        y_t = b0 + phi * y_t-1 - phi*b0 + b1*light_t + w_t
#                   y_t = b0*(1 - phi) + phi * y_t-1 + b1*light_t + w_t
# now with the centered data, the intercept is defined in terms of phi
for (i in 1:365){
  z[i+1] = b0 * (1 - phi) + z[i] * phi + light.rel[i] * b1 + rnorm(1, 0, sdp)
}

# Estimate observed data
# set.seed(4)
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
sink("ss_light_centered.stan")

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
       z[i] ~ normal(b0*(1-phi)+ z[i-1]*phi+light[i]*b1, sdp);// process model with error
        }


 }

      "

    ,fill=TRUE)
sink()
closeAllConnections()

#Prep model
model<-"ss_light_centered.stan"
model<-stan_model(model)

#Create data object
data <- list(N = length(y[1:365]),y=y[1:365], z0=y[1], light=light.rel)

#Run Stan
fit_cent<-rstan::sampling(object=model, data = data,  iter = 3000, chains = 4)#, control=list(max_treedepth = 15))


##Print and extract
fit_extract<-rstan::extract(fit_cent)
print(fit_cent, pars=c( "sdo","phi", "sdp", "b0","b1" ))

pairs(fit_cent, pars=c("sdo","phi", "sdp","b0","b1","lp__"))

# This definitely helps with the b0/phi equifinality, while not as much with
# the b1/phi and b1/b0, although it seems like both of these relationships
# are not as strong as they were in the uncentered model.

##HMC diagnostics
rstan::check_hmc_diagnostics(fit_cent)

##Traceplots

traceplot(fit_cent, pars=c("sdo","phi", "sdp", "b1", "b0"))

##Plot density plots
plot_sdo <- stan_dens(fit_cent, pars="sdo") + geom_vline(xintercept =sdo)
plot_sdp <- stan_dens(fit_cent, pars="sdp") + geom_vline(xintercept = sdp)
plot_phi <- stan_dens(fit_cent, pars="phi") + geom_vline(xintercept = phi)
plot_b0 <- stan_dens(fit_cent, pars="b0") + geom_vline(xintercept = 2)+xlab("Intercept")
plot_b1 <- stan_dens(fit_cent, pars="b1") + geom_vline(xintercept = b1)+xlab("Light beta")
grid.arrange(plot_b0,plot_b1,plot_phi,plot_sdp,plot_sdo, nrow=2)

##Create object with estimated missing data

fit_summary<-summary(fit_cent, probs=c(0.05,.5,.95))$summary
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






