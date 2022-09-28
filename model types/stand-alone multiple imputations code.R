
# load libraries ----------------------------------------------------------

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(data.table)
library(Amelia)
library(plyr)


# set wd ------------------------------------------------------------------

setwd("~/GitHub/missing-data")


# Simulate data ----------------------------------------------------------

#basic parameters
N<-365 #length of data
z<-numeric(N+1) 
ts<-seq(1:N)
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
#Run the function
light<-lightest(time, 47.8762, -114.03, 105) #location is Flathead Bio Station just for fun
light.rel<-light/max(light) #make relative light



#Simulate GPP based on time-series model with known parameters
#set parameters
z<-NA
sdp <- 0.1 # Process error. Must be >0 and <1
sdo<-0.1 # Observation error. Must be >0 and <1. Problems arise when sdo > sdp
phi <-0.8 #Auto-correlation. Must be >0 and <1
b0<-0 # Intercept. No limits, but generally positive. Acts as scaling factor for simulated data (changes the magnitude) 
b1<-0.6 #Light coefficient. Defines magnitude of hump shape.
z[1]<-1 # Initial value of state. May need to change to make initial simulations smoother. 

## Set the seed, so we can reproduce the results
set.seed(5)
## For-loop that simulates the state through time.
for (i in 1:N){
  z[i+1] = phi*z[i]+light.rel[i]*b1+rnorm(1, 0, sdp)
}

#Estimate observed data
set.seed(4)
sd.o <-rnorm(N, 0, sdo)
y <-z[2:(N+1)]+ sd.o
sdo_sd<-mean(abs(sdo))

##Bind simulated data and time
y_full<-as.data.frame(cbind(z[1:365],y,time))
names(y_full)<-c("z", "y", "time")

ggplot(data=y_full, aes(x=time, y=z))+
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




# Force some missing data ------------------------------------------
#Note-the model will still work if no data is missing but the objects in the next block of code (i.e., "y_index_mis", etc.)
#still need to be created

##Set seed for replicability
set.seed(50)

##Load data from SDW characterizing lengths
obs<-readRDS("C:/Users/matt/IDrive-Sync/Postdoc/Estimating missing data/daily_predictions/SDW_obs_length.RDS")
miss<-readRDS("C:/Users/matt/IDrive-Sync/Postdoc/Estimating missing data/daily_predictions/SDW_miss_length.RDS")

##Replicate simulated data
n.datasets<-41
y_miss<-list(y)
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
miss.num<-sapply(y_miss, function(y) sum(length(which(is.na(y)))))
prop.miss<-round(miss.num/length(y)*100)
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


# Use Amelia to impute missing data ---------------------------------------

#single data set

#create data set with missing data
y_miss1<-na_if(y_miss[[3]], -100)
w<-as.data.frame(cbind(y_miss1,light,ts))# Combine GPP, light, and ts to one data frame for Amelia

# run amelia
data.amelia<-amelia(w, m = 5, p2s=1, ts="ts", lags="y_miss1") #m=number of imputations. ts= time component, lags= which column has a 1-day autocorrelation, 
summary(data.amelia)


#Plot imputed and 'observed' (simulated) data
  z1<-data.amelia[[1]]$imp1$y_miss1
  z2<-data.amelia[[1]]$imp2$y_miss1
  z3<-data.amelia[[1]]$imp3$y_miss1
  z4<-data.amelia[[1]]$imp4$y_miss1
  z5<-data.amelia[[1]]$imp5$y_miss1
  q1<-z1[y_index_mis[[3]]]
  q2<-z2[y_index_mis[[3]]]
  q3<-z3[y_index_mis[[3]]]
  q4<-z4[y_index_mis[[3]]]
  q5<-z5[y_index_mis[[3]]]
  q.ts<-ts[y_index_mis[[3]]]
  q<-cbind(q1,q2,q3,q4,q5)
  q.mean<- apply(q[,-1], 1, mean)
  q.min<-apply(q[,-1], 1, min)
  q.max<-apply(q[,-1], 1, max)
  q.sd<- apply(q[,-1], 1, sd)
  x_obs_data<-as.data.frame(cbind(y_miss[[1]],ts))
  x_imp_data<-as_tibble(cbind(q.mean,q.sd,q.ts, q.min,q.max))
  x_imp_data<-as.data.frame(x_imp_data %>% slice(1:length(y_index_mis[[3]])))
  x_obs<-y_miss[[1]][y_index_mis[[3]]]
  x_combined<-as.data.frame(cbind(x_imp_data$q.mean,x_imp_data$q.sd,x_obs,x_imp_data$q.min,x_imp_data$q.max))
  ggplot(data=x_imp_data, aes(x=q.ts, y=q.mean))+
    geom_errorbar(data=x_imp_data,aes(x=q.ts,ymin=q.min, ymax=q.max))+
    geom_point(data=x_obs_data, aes(x=ts, y=V1, color="gray"), size=3)+
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
  
  ggplot(data=x_combined, aes(x=x_obs, y=V1))+
    geom_abline(intercept=0, slope=1)+
    theme_classic()+
    geom_errorbar(aes(ymin=q.min, ymax=q.max))+
    geom_point( size=3)+
    theme(legend.position="top")+
    ylab("Imputed data (mean and range)")+
    xlab("Simulated data")+
    theme(axis.title.x=element_text(size=18,colour = "black"))+
    theme(axis.title.y=element_text(size=18,colour = "black"))+
    theme(axis.text.y=element_text(size=18,colour = "black"))+
    theme(axis.text.x=element_text(size=18,colour = "black"))
  compare.density(data.amelia, var = "y_miss1", xlab="Data", main="", lwd=2)
  
  # STAN code ---------------------------------------------------------------
  
  
  sink("miss_ss.stan")
  
  cat("
 
/*----------------------- Data --------------------------*/
  /* Data block: defines the objects that will be inputted as data */
  data {
    int N; // Length of state and observation time series
    vector[N] y_miss; // Observations
    real z0; // Initial state value
    vector[N] light;
    vector[N] miss_vec; // index or location of missing values within the dataset
  }
  
/*----------------------- Parameters --------------------------*/
  /* Parameter block: defines the variables that will be sampled */
  parameters {
    real<lower=0> sdp; // Standard deviation of the process equation
    real<lower=0> sdo; // Standard deviation of the observation equation
    real b0;
    real b1;
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    vector[N] z; // State time series
  }
  

  /*----------------------- Model --------------------------*/
  /* Model block: defines the model */
  model {
  int count;
    // Prior distributions
    sdo ~ normal(0.1, 0.01);
    sdp ~ normal(0, 1);
    phi ~ beta(1,1);
    b0 ~ normal(0,5);
    b1 ~ normal(0,5);
    
    // Distribution for the first state
   z[1] ~ normal(z0, sdp);
   
   
    // Distributions for all other states
    count = 0;
    for(t in 2:N){
       z[t] ~ normal(b0+ z[t-1]*phi+light[t]*b1, sdp);
    if(miss_vec[t]==0){
      y_miss[t] ~ normal(z[t], sdo); // observation model with fixed observation error
    }
    }
  
  }
  
    "
      ,fill=TRUE)
  sink()
  
  
  # MI-model (with loop) ----------------------------------------------------
  
  
  fit.amelia <- vector("list",1)
  data.stan.amelia <- vector("list",(length(1)))
  list.1<- vector("list",(5))
  for(i in 1){##two because the list holding the lists starts with a list of misc data that arent the estimates
    for(g in 1:5){
      list.1[[g]] <- data.amelia[["imputations"]][[g]]$y_miss1
      
    }
    data.stan.amelia[[i]]<-list.1
  }
  
  plot(data.stan.amelia[[1]][[2]],data.stan.amelia[[1]][[3]])
  
  ##Run Stan
  
  
  fit.stan.miss.amelia <- vector("list",1)
  fit.stan.miss.imp<- vector("list",5)
  model<-"miss_ss.stan"
  model<-stan_model(model)
  for(i in 3){
    for(g in 1:5){
      ##Load data
      data <- list(  N = length(y_miss[[i]]),
                     y_nMiss = unlist(y_nMiss[[1]]),
                     y_index_mis = unlist(y_index_mis[[1]]),
                     y_miss= unlist(data.stan.amelia[[1]][[g]]),
                     light=light.rel,
                     z0=y_miss[[i]][1],
                     miss_vec=miss_vec[[3]])
      
      ##Run Stan
      fit.stan.miss.imp[[g]]<- rstan::sampling(object=model, data = data,  iter = 4000, chains = 4, control=list(max_treedepth = 15))
      print(g)
      rstan::check_hmc_diagnostics(fit.stan.miss.imp[[g]])
      }
    fit.stan.miss.amelia[[i]]<-fit.stan.miss.imp
  
  }
  
  # MI-write or read RDS files ----------------------------------------------
  
  #saveRDS(fit.stan.miss.amelia, file = "test.RDS")
  fit.stan.miss.amelia<-readRDS("test.RDS")
  
  # MI-pull and clean paramters ---------------------------------------------
  
  
  ##Pull param estimates into list
  fit_summary_pars_amelia <- vector("list",1)
  list.2<- vector("list",1)
  for (i in 3){
    for (g in 1:5){
      list.2[[g]] <-summary(fit.stan.miss.amelia[[i]][[g]], pars=c("sdp","sdo","phi", "b1", "b0"), probs=c(0.025,.5,.975))$summary
    }
    fit_summary_pars_amelia[[i]]<-list.2
  }
  fit_summary_pars_amelia[[3]]
  
  ##Unlist,cleanup, and add factors
  fit_summary<- NA
  
    fit_summary<-as.data.frame(do.call(rbind,fit_summary_pars_amelia[[3]]))
    fit_summary$param<-rep(c("sdp","sdo","phi", "b1","b0"), times=length(1) )#add parameter factor
    fit_summary$prop.missing<-rep(prop.miss[3], each=1) #add prop of missing data
    fit_summary$m<-rep(c(1:5), times=5)
    fit_summary$m<-as.factor(fit_summary$m)
    row.names(fit_summary)<-1:25 #remove row names
    fit_summary$param<-as.factor(fit_summary$param) #make factor
    fit_summary$prop.missing<-as.factor(fit_summary$prop.missing) #make factor
    
    colnames(fit_summary)<-c("mean", "se_mean", "sd", "min", "med", "high","n_eff", "rhat", "param","prop.missing", "m")

  
##Check that it looks good
  summary(fit_summary) #check it looks good
  head(fit_summary)

##Model average based on Rubins Rules
#Calculate the average and within imputation variation
summ<-ddply(fit_summary, c("param", "prop.missing"), summarize,
              avg= mean(mean),
              vw=sum(sd^2)/5)

#join to full summary for further calcs
fit_summary<-left_join(fit_summary, summ)
  
#calculate square difference of mean and individual estimates
summ2<-ddply(fit_summary, c("param","m"), summarize,
             square.diff= (mean-avg)^2)

#Add that to main summary
fit_summary$square.diff<-summ2$square.diff

#calculate between imputation variance
summ3<-ddply(fit_summary, c("param"), summarize,
            vb= sum(square.diff)/4)

#Add to main summary object
summ$vb<-summ3$vb

#calculate total variance
summ$vt<-summ$vw+summ$vb+(summ$vb/5)

#calculate standard error
summ$se<-sqrt(summ$vt)

#PLOT
#create known data object for plotting
known.data<-c(sdp,sdo,phi,b1,b0)
  known<-as.data.frame(known.data)
  
  known.param<-c("sdp","sdo","phi", "b1","b0")
  #known.missing<-c(5, 5, 5,5)
  known<-as.data.frame(cbind(known.data, known.param))
  known$known.data<-as.numeric(as.character(known$known.data))
  
  ggplot(data=summ, aes(x=avg, y=prop.missing ))+
    geom_point(aes( color=param, group=param),size=3,position=position_dodge(0.5))+
    #geom_point(data=fit_summary[[3]],aes(x=mean, y=prop.missing, color=param, group=param),size=3,position=position_dodge(0.5))+
    #geom_point(data=fit_summary[[4]],aes(x=mean, y=prop.missing, color=param, group=param),size=3,position=position_dodge(0.5))+
    theme_classic()+
    geom_errorbar(aes(xmin=avg-se, xmax=avg+se,color=param, group=param), width=0.2,position=position_dodge(0.5))+
    theme(legend.position="top")+
    ylab("Percent of Missing Data")+
    xlab("Parameter Estimate")+
    scale_color_manual(values=c("blue", "green", "darkgray", "black", "red"))+
    theme(axis.title.x=element_text(size=18,colour = "black"))+
    theme(axis.title.y=element_text(size=18,colour = "black"))+
    theme(axis.text.y=element_text(size=18,colour = "black"))+
    theme(axis.text.x=element_text(size=18,colour = "black"))+
    geom_vline(xintercept = known.data,color=c("red", "black", "darkgray", "green","blue"))
  
  
  
  
#MANY data sets

#Impute
data.amelia<- vector("list",max(n.datasets)) #Create list of datasets
for(i in 2:max(n.datasets)){ #start at 2 becaue the first dataset does not have missing data
  y_miss1<-na_if(y_miss[[i]], -100)
  w<-as.data.frame(cbind(y_miss1,light,ts))# Combine GPP, light, and ts to one data frame for Amelia
  z2<-amelia(w, m = 5, p2s=1, ts="ts", lags="y_miss1") #m=number of imputations. ts= time component, lags= which column has a 1-day autocorrelation, 
  summary(z2)
  data.amelia[[i]]<-z2
}


#Assign NAs to the missing data in each list
for(i in 2:41){
  z1<-data.amelia[[i]]$imputations$imp1$y_miss1
  z2<-data.amelia[[i]]$imputations$imp2$y_miss1
  z3<-data.amelia[[i]]$imputations$imp3$y_miss1
  z4<-data.amelia[[i]]$imputations$imp4$y_miss1
  z5<-data.amelia[[i]]$imputations$imp5$y_miss1
  q1<-z1[y_index_mis[[i]]]
  q2<-z2[y_index_mis[[i]]]
  q3<-z3[y_index_mis[[i]]]
  q4<-z4[y_index_mis[[i]]]
  q5<-z5[y_index_mis[[i]]]
  q.ts<-ts[y_index_mis[[i]]]
  q<-cbind(q1,q2,q3,q4,q5)
  q.mean<- apply(q[,-1], 1, mean)
  q.min<-apply(q[,-1], 1, min)
  q.max<-apply(q[,-1], 1, max)
  q.sd<- apply(q[,-1], 1, sd)
  x_obs_data<-as.data.frame(cbind(y_full$y,ts))
  x_imp_data<-as_tibble(cbind(q.mean,q.sd,q.ts, q.min,q.max))
  x_imp_data<-as.data.frame(x_imp_data %>% slice(1:length(y_index_mis[[i]])))
  x_obs<-y_full$y[y_index_mis[[i]]]
  x_combined<-as.data.frame(cbind(x_imp_data$q.mean,x_imp_data$q.sd,x_obs,x_imp_data$q.min,x_imp_data$q.max))
 g<-ggplot(data=x_imp_data, aes(x=q.ts, y=q.mean))+
    geom_errorbar(data=x_imp_data,aes(x=q.ts,ymin=q.min, ymax=q.max))+
    geom_point(data=x_obs_data, aes(x=ts, y=V1, color="gray"), size=3)+
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
  
  s<-ggplot(data=x_combined, aes(x=x_obs, y=V1))+
    geom_abline(intercept=0, slope=1)+
    theme_classic()+
    geom_errorbar(aes(ymin=q.min, ymax=q.max))+
    geom_point( size=3)+
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



# MI-write or read RDS files ----------------------------------------------

saveRDS(fit.stan.miss.amelia, file = "full_sim_week_amelia_sdp_01_phi_8_b0_1_b1_1.RDS")
saveRDS(fit_summary, file = "summary_sim_week_amelia_sdp_01_phi_8_b0_1_b1_1.RDS")
saveRDS(known, file = "known_sim_week_amelia_sdp_01_phi_8_b0_1_b1_1.RDS")
