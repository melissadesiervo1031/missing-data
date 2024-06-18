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
library(Amelia)
library(plyr)

# Simulate data -----------------------------------------------------------

N<-365 #length of data
t<-0:(N-1) #time
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
light.c<-(light.l-mean(light.l))
mean(light.c)
sd(light.c)


##GPP based on time-series model with known paramters
set.seed(553)
x<-NA
sdp <- 0.1
phi <-0.8
b0<-0.1
b1<-0.1
x[1]<-3.5

##GPP based on time-series model with known paramters
set.seed(553)
x<-NA
sdp <- 0.01
phi <-0.2
b0<-0.1
b1<-0.1
x[1]<-0.9

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
y_miss_40<-as.data.frame(cbind(y_miss[[8]],t))

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
  ggtitle("40% removed")



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


# Force some missing data-BY DAY ------------------------------------------

#Note-the model will still work if no data is missing but the objects in the next block of code (i.e., "y_index_mis", etc.)
#still need to be created
set.seed(90)
#Center
y_miss<-list(x,x,x,x,x,x,x,x)

#vector of missing number amounts
missing_n<-c(0,round(N*0.15),round(N*0.3),round(N*0.4),round(N*0.5),round(N*0.6),round(N*0.7),round(N*0.8))
#Create a list to store missing data integers
y_missing_integers<-c()
#Create a vector of the length of data representing the indexes
index<-c()
index[[1]]<-1:N

#Create initial missing data integer set. Then build on with the for loop
y_missing_integers[[1]]<-sample(index,missing_n[1],replace = FALSE)
#index<-index[-y_missing_integers[[1]]]

#Adds the previous missing data integers to the next set. This makes them nested.
#i.e, all the missing data 10 integers are added to 15 more to make missing data 15.
for(i in 2:length(missing_n)){
 h<-sample(index[[i-1]],(missing_n[i]-missing_n[i-1]), replace = FALSE)
 index[[i]]<-setdiff(index[[i-1]],h)
 y_missing_integers[[i]]<-unlist(c(h,y_missing_integers[[i-1]]))
}

#Assign NAs to the missing data in each list
for(i in 2:length(missing_n)){
  z<-y_miss[[i]]
  z[y_missing_integers[[i]]]<-NA
  y_miss[[i]]<-z
}

#y_miss.o<-y_miss
#Center
#for(i in 1:length(missing_n)){
 # m<-mean(y_miss[[i]],na.rm=TRUE) 
#  s<-sd(y_miss[[i]],na.rm=TRUE)
 # y_miss[[i]]<-(y_miss[[i]]-m)/s
#}

#Check that each list has the right amount of missing data (or close to it...)
miss<-sapply(y_miss, function(x) sum(length(which(is.na(x)))))
prop.miss<-round(miss/length(x)*100)
prop.miss

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
for(i in 2:length(missing_n)){
  r<-y_miss[[i]]
  r[y_missing_integers[[i]]]<--100 #arbitrary number
  r[1]<-x[1]
  #r[1]<-(x[1]-mean(y_miss.o[[i]], na.rm=TRUE))/sd(y_miss.o[[i]], na.rm=TRUE)
  y_miss[[i]]<-r
} 

summary(y_miss)



# Stan Code ---------------------------------------------------------------


sink("toy2p.stan")

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

sink("toy2p_cent.stan")
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
  transformed data {
  vector[N] light_std;
  vector[N] y_std;
  light_std = (light - mean(light)) / sd(light);
  y_std = (y_miss - mean(y_miss)) / sd(y_miss);
}
  
/*----------------------- Parameters --------------------------*/
  /* Parameter block: defines the variables that will be sampled */
  parameters {
    vector<lower = 0>[y_nMiss] y_imp;// Missing data
    real<lower=0> sdp_std; // Standard deviation of the process equation
    real b0_std;
    real b1_std;
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    
  }
  transformed parameters { 
    vector[N] y;
    y=y_std; // makes the data a transformed variable
    y[y_index_mis] =y_imp; // replaces missing data in y with estimated data
    } 
  /*----------------------- Model --------------------------*/
  /* Model block: defines the model */
  model {
    // Prior distributions
   
    sdp_std ~ normal(0, 1);
    phi ~ beta(1,1);
    b0_std ~ normal(0,5);
    b1_std ~ normal(0,5);
    
    // Distribution for the first state
    y_std[1] ~ normal(z0, sdp_std);
    
    // Distributions for all other states
    for(t in 2:N){
      y_std[t] ~ normal(b0_std+y_std[t-1]*phi+light_std[t]*b1_std, sdp_std);
    }
    
  }
  generated quantities {
 vector[N] y_rep; // replications from posterior predictive dist
 real b0;
 real b1;
 real<lower=0> sdp;
 b0 = sd(y_miss) * (b0_std - (b1_std * mean(light)/sd(light)))+mean(y_miss);
 b1 = b1_std * sd(y_miss)/ sd(light);
 sdp= sd(y_miss)* sdp_std;
 
 y_rep[1]=normal_rng(y[1], 0.1);
 
 for (t in 2:N) {
 y_rep[t]=normal_rng(b0+y_miss[t-1]*phi+light[t]*b1, sdp);
 
 }
 
  }

    "
,fill=TRUE)
sink()
closeAllConnections()

# Bayes parameter estimation model ----------------------------------------

fit.miss <- vector("list",length(missing_n))
model<-"toy2p.stan"
model<-stan_model(model)

#for(i in 1:8){
  for(i in 1){
##Load data
data <- list(   N = length(y_miss[[i]]),
                y_nMiss = unlist(y_nMiss[[i]]),
                y_index_mis = unlist(y_index_mis[[i]]),
                y_miss=unlist(y_miss[[i]]) ,
                light=(light.l-mean(light.l))/sd(light.l),
                z0=y_miss[[i]][1])

##Run Stan
fit.stan.miss.amelia[[i]]<- rstan::sampling(object=model, data = data,  iter = 4000, chains = 4)#, control=list(max_treedepth = 15,adapt_delta = 0.99))
  }


##Increase step size
#adapt_delta = 0.99
print(fit.miss[[2]])

##Pull param estimates into list
fit_summary_pars_bayes <- vector("list",length(missing_n))
for (i in 1){
  fit.stan.miss.amelia[[i]]<-(summary(fit.stan.miss.amelia[[i]], pars=c("sdp","phi", "b1", "b0"), probs=c(0.025,.5,.975))$summary)
   }
fit_summary_pars_bayes[[1]]


fit_summary_bayes<- vector("list",length(missing_n))
##Unlist,cleanup, and add factors
fit_summary_bayes<-as.data.frame(do.call(rbind.data.frame, fit_summary_pars_bayes)) #Unlist
fit_summary_bayes$param<-rep(c("sdp","phi", "b1","b0"), times=length(missing_n) )#add parameter factor
fit_summary_bayes$prop.missing<-rep(prop.miss, times=c(4,4,4,4,4,4,4,4)) #add prop of missing data
row.names(fit_summary_bayes)<-1:(length(missing_n)*4) #remove row names
fit_summary_bayes$param<-as.factor(fit_summary_bayes$param) #make factor
fit_summary_bayes$prop.missing<-as.factor(fit_summary_bayes$prop.missing) #make factor
colnames(fit_summary_bayes)<-c("mean", "se-mean", "sd", "min", "med", "high","n_eff", "rhat", "param","prop.missing")
summary(fit_summary_bayes) #check it looks good
head(fit_summary_bayes)
fit_summary_bayes<-as.data.frame(fit_summary_bayes)

known.data<-c(sdp,phi,b1,b0)
known<-as.data.frame(known.data)

known.param<-c("sdp","phi", "b1","b0")
#known.missing<-c(5, 5, 5,5)
known<-as.data.frame(cbind(known.data, known.param))
known$known.data<-as.numeric(as.character(known$known.data))

ggplot(data=fit_summary_bayes, aes(x=mean, y=prop.missing ))+
  geom_point(aes( color=param, group=param),size=3,position=position_dodge(0.5))+
  theme_classic()+
  geom_errorbar(aes(xmin=min, xmax=high,color=param, group=param), width=0.2, size=0.5,position=position_dodge(0.5))+
  theme(legend.position="top")+
  ylab("Percent of Missing Data")+
  xlab("Parameter Estimate")+
  scale_color_manual(values=c("blue", "green", "darkgray", "black"))+
  scale_x_continuous(limits=c(-90,60))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  geom_vline(xintercept = known.data,color=c("black","darkgray","green","blue"))
  

saveRDS(fit.miss, file = "full_sim_day_bayes_sdp_01_phi_2_b0_1_b1_1.RDS")
saveRDS(fit_summary_bayes, file = "summary_sim_day_bayes_sdp_01_phi_2_b0_1_b1_1.RDS")
saveRDS(known, file = "known_sim_day_bayes_sdp_01_phi_2_b0_1_b1_1.RDS")


write.csv(fit_summary_pars_bayes, "Rstan_missing data_bayes_sdp_01_phi_5.csv")
##Print and extract
fit_extract<-rstan::extract(fit[[1]])
print(fit.miss, pars=c("sdp","phi", "b1","b0"))

pairs(fit.miss[[5]], pars=c("sdp","phi", "b1", "b0"))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit.miss[[2]])

##Traceplots

traceplot(fit.miss[[2]], pars=c("phi", "b0","sdp", "b1"))

##Posterior densities compared to known parameters
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

fit_summary<-summary(fit.miss[[2]], probs=c(0.025,.5,.975))$summary

my_data <- as_tibble(fit_summary)
y_est_data<-my_data %>% slice(1:length(y_index_mis[[2]]))

##Create object with observed data
date<-1:N
y_obs_data<-as.data.frame(cbind(date[y_index_mis[[2]]],x[y_index_mis[[2]]]))
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
  geom_errorbar(aes(ymin=low, ymax=high))+
  theme(legend.position="top")+
  ylab("Estimated GPP")+
  xlab("Observed GPP")+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))


raw_data<-as.data.frame(cbind(x,date))

ggplot(data=y_combined, aes(x=date, y=estimated))+
  geom_point(data=raw_data, aes(x=date, y=x, color="gray"), size=3)+
  geom_errorbar(data=y_combined,aes(x=date,ymin=low, ymax=high))+ 
  geom_point(aes(color="red"), size=3)+
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  scale_color_identity(guide = "legend", name=element_blank(), labels=c("Simulated data", "Imputed data"))+
  theme(legend.position="top")+
  ylab("Data")+
  xlab("Time")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



# Missing Imputations model -----------------------------------------------

ts<-as.Date(time,format="%m/%d/%Y")

data.amelia<- vector("list",length(missing_n))
for(i in 1:length(missing_n)){
y_miss1<-na_if(y_miss[[i]], -100)
y_miss1[i]<-x[1]
#y_miss1[1]<-(x[1]-mean(y_miss.o[[i]], na.rm=TRUE))/sd(y_miss.o[[i]], na.rm=TRUE)
w<-as.data.frame(cbind(y_miss1,light,ts))
z2<-amelia(w, m = 3, p2s=1, ts="ts", lags="y_miss1")
summary(z2)
data.amelia[[i]]<-z2
}


#Assign NAs to the missing data in each list
for(i in 6){
  z1<-data.amelia[[i]]$imputations$imp1$y_miss1
  z2<-data.amelia[[i]]$imputations$imp2$y_miss1
  z3<-data.amelia[[i]]$imputations$imp3$y_miss1
  z4<-data.amelia[[i]]$imputations$imp4$y_miss1
  z5<-data.amelia[[i]]$imputations$imp5$y_miss1
  q1<-z1[y_missing_integers[[i]]]
  q2<-z2[y_missing_integers[[i]]]
  q3<-z3[y_missing_integers[[i]]]
  q4<-z4[y_missing_integers[[i]]]
  q5<-z5[y_missing_integers[[i]]]
  q.ts<-ts[y_missing_integers[[i]]]
  q<-cbind(q1,q2,q3,q4,q5)
  q.mean<- apply(q[,-1], 1, mean)
  q.min<-apply(q[,-1], 1, min)
  q.max<-apply(q[,-1], 1, max)
  q.sd<- apply(q[,-1], 1, sd)
  x_obs_data<-as.data.frame(cbind(x,ts))
  x_imp_data<-as_tibble(cbind(q.mean,q.sd,q.ts, q.min,q.max))
  x_imp_data<-as.data.frame(x_imp_data %>% slice(1:length(y_missing_integers[[i]])))
  x_obs<-x[y_missing_integers[[i]]]
  x_combined<-as.data.frame(cbind(x_imp_data$q.mean,x_imp_data$q.sd,x_obs,x_imp_data$q.min,x_imp_data$q.max))
  g<-ggplot(data=x_imp_data, aes(x=q.ts, y=q.mean))+
    geom_point(data=x_obs_data, aes(x=ts, y=x, color="gray"), size=3)+
    geom_point(aes(color="red"), size=3)+
    geom_errorbar(data=x_imp_data,aes(x=q.ts,ymin=q.min, ymax=q.max))+ 
    theme_classic()+
    theme(axis.title.x=element_text(size=18,colour = "black"))+
    theme(axis.title.y=element_text(size=18,colour = "black"))+
    theme(axis.text.y=element_text(size=18,colour = "black"))+
    theme(axis.text.x=element_text(size=18,colour = "black"))+
    scale_color_identity(guide = "legend", name=element_blank(), labels=c("Simulated data", "Imputed data"))+
    theme(legend.position="top")+
    ylab("Data")+
    xlab("Time")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  s<-ggplot(data=x_combined, aes(x=x_obs, y=V1))+
    geom_point( size=3)+
    geom_abline(intercept=0, slope=1)+
    theme_classic()+
    geom_errorbar(aes(ymin=q.min, ymax=q.max))+
    theme(legend.position="top")+
    ylab("Imputed data (mean and range)")+
    xlab("Simulated data")+
    theme(axis.title.x=element_text(size=18,colour = "black"))+
    theme(axis.title.y=element_text(size=18,colour = "black"))+
    theme(axis.text.y=element_text(size=18,colour = "black"))+
    theme(axis.text.x=element_text(size=18,colour = "black"))
  compare.density(data.amelia[[i]], var = "y_miss1", xlab="Data", main="", lwd=2)
  plot(g)
  plot(s)
  
}
     

?compare.density


fit.amelia <- vector("list",length(missing_n))
data.stan.amelia <- vector("list",(length(missing_n)))
list.1<- vector("list",(3))
for(i in 1:length(missing_n)){##two because the list holding the lists starts with a list of misc data that arent the estimates
  for(g in 1:3){
 list.1[[g]] <- data.amelia[[i]]$imputations[[g]]$y_miss1
         
}
  data.stan.amelia[[i]]<-list.1
  }

plot(data.stan.amelia[[8]][[2]],data.stan.amelia[[8]][[3]])
##Run Stan


  fit.stan.miss.amelia <- vector("list",length(missing_n))
  fit.stan.miss.imp<- vector("list",3)
  model<-"toy2p_cent.stan"
  model<-stan_model(model)
  for(i in 2:8){
    for(g in 1:3){
    ##Load data
    data <- list(  N = length(y_miss[[1]]),
                   y_nMiss = unlist(y_nMiss[[1]]),
                   y_index_mis = unlist(y_index_mis[[1]]),
                   y_miss= unlist(data.stan.amelia[[i]][[g]]),
                   light=light.l,
                   z0=unlist(data.stan.amelia[[i]][[g]])[1])
    
    ##Run Stan
    fit.stan.miss.imp[[g]]<- rstan::sampling(object=model, data = data,  iter = 4000, chains = 4, control=list(max_treedepth = 15))
  }
    fit.stan.miss.amelia[[i]]<-fit.stan.miss.imp
  }
    print(fit.stan.miss.amelia[[2]][[1]])
  ##Pull param estimates into list
  fit_summary_pars_amelia <- vector("list",length(missing_n))
  list.2<- vector("list",3)
  
  for (i in 2:8){#:length(missing_n_week)){
    for (g in 1:3){
      list.2[[g]] <-summary(fit.stan.miss.amelia[[i]][[g]], pars=c("sdp","phi", "b1", "b0"), probs=c(0.025,.5,.975))$summary
  }
  fit_summary_pars_amelia[[i]]<-list.2
  }
  
  fit_summary_pars_amelia[[8]][[1]]
  
  ##Unlist,cleanup, and add factors
  fit_summary<- vector("list",3)
  for(i in 2:8){
    fit_summary[[i]]<-as.data.frame(do.call(rbind,fit_summary_pars_amelia[[i]]))
  fit_summary[[i]]$param<-rep(c("sdp","phi", "b1","b0"), times=length(3) )#add parameter factor
  fit_summary[[i]]$prop.missing<-rep(prop.miss[i], each=4,4,4,4,4,4,4,4) #add prop of missing data
  row.names(fit_summary[[i]])<-1:12 #remove row names
  fit_summary[[i]]$param<-as.factor(fit_summary[[i]]$param) #make factor
  fit_summary[[i]]$prop.missing<-as.factor(fit_summary[[i]]$prop.missing) #make factor
  colnames(fit_summary[[i]])<-c("mean", "se.mean", "sd", "min", "med", "high","n_eff", "rhat", "param","prop.missing")
  }
  fit_summary<-as.data.frame(do.call(rbind,fit_summary))
  fit_summary_avg <- ddply(fit_summary, c("param", "prop.missing"), summarize,
                 average = mean(mean))
                 
  fit_summary_avg<-as.data.frame(fit_summary_avg)
  
  summary(fit_summary) #check it looks good
 head(fit_summary)
  known.data<-c(sdp,phi,b1,b0)
  known<-as.data.frame(known.data)
  
  known.param<-c("sdp","phi", "b1","b0")
  #known.missing<-c(5, 5, 5,5)
  known<-as.data.frame(cbind(known.data, known.param))
  known$known.data<-as.numeric(as.character(known$known.data))
  
  ggplot(data=fit_summary, aes(x=mean, y=prop.missing ))+
    geom_point(aes( color=param, group=param),size=3,position=position_dodge(0.5))+
    theme_classic()+
    geom_errorbar(data=fit_summary, aes(xmin=min, xmax=high,color=param, group=param), width=0.2,position=position_dodge(0.5))+
    geom_point(data=fit_summary_avg,aes(x=average,y=prop.missing, color=param, group=param),size=3,position=position_dodge(0.5),shape=2)+
    theme(legend.position="top")+
    ylab("Percent of Missing Data")+
    xlab("Parameter Estimate")+
    scale_color_manual(values=c("blue", "green", "darkgray", "black"))+
    scale_x_continuous(limits=c(-.1,0.9))+
    theme(axis.title.x=element_text(size=18,colour = "black"))+
    theme(axis.title.y=element_text(size=18,colour = "black"))+
    theme(axis.text.y=element_text(size=18,colour = "black"))+
    theme(axis.text.x=element_text(size=18,colour = "black"))+
    geom_vline(xintercept = known.data,color=c("blue", "darkgray", "green", "black"))
  
  #saveRDS(fit.stan.miss.amelia, file = "full_sim_day_amelia.RDS")
  #saveRDS(fit_summary, file = "summary_sim_day_amelia.RDS")
  #saveRDS(known, file = "known_sim_day_amelia.RDS")
  