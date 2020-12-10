
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(dplyr)
library(shinystan)
library(faux)
###Simulate data

##Simulate nitrogen, srp, GPP, and light with a degree of correlation and random noise (process error) centered at 0
# Set the seed, so we can reproduce the results
set.seed(620)
dat1 <- rnorm_multi(n = 366,
                   mu = c(0,0,0,0),
                   sd = c(5,5,0.5,0.5),
                   r = c(1,.6,.1,.3,.6,1,.5,.5,.1,.5,1,.9,.3,.5,.9,1), ##correlation matrix
                   varnames = c("srp", "nit", "light", "GPP"),
                   empirical = FALSE)
dat1$Y<-rep(1,366)
quad<-abs(seq(-2,2, length=366))
dat1$light<-dat1$light-quad
plot(time, dat1$light)

dat2 <- rnorm_multi(n = 366,
                    mu = c(0,0,0,0),
                    sd = c(5,5,0.5,0.5),
                    r = c(1,.4,.1,.3,.4,1,.5,.8,.1,.5,1,.9,.3,.8,.9,1), ##correlation matrix
                    varnames = c("srp", "nit", "light", "GPP"),
                    empirical = FALSE)
dat2$Y<-rep(2,366)
quad<-abs(seq(-2,2, length=366))
dat2$light<-dat2$light-quad
plot(time, dat2$light)


dat3 <- rnorm_multi(n = 366,
                    mu = c(0,0,0,0),
                    sd = c(5,5,0.5,0.5),
                    r = c(1,.4,.1,.3,.4,1,.5,.8,.1,.5,1,.6,.3,.8,.6,1), ##correlation matrix
                    varnames = c("srp", "nit", "light", "GPP"),
                    empirical = FALSE)

dat3$Y<-rep(3,366)
quad<-abs(seq(-5,5, length=366))
dat3$light<-dat3$light-quad
timedat3<-1:366
plot(timedat3, dat3$light)

data<-rbind(dat1,dat2,dat3)
nit<-dat$nit
srp<-dat$srp


###Create observed data
# Create a vector that will keep track of the observations
# It's of length T
y <- numeric(TT)
# Standard deviation of the observation error
sdo <- 1
# For t=1, ... T, add measurement error
# Remember that z[1] is t=0
y <- data$GPP[2:(TT+1)] + rnorm(TT, 0, sdo)



###Plot simulated and observed data
plot(1:TT, y,
     pch=3, cex = 0.8, col="blue", ty="o", lty = 3,
     xlab = "t", ylab = expression(z[t]),
     #xlim = c(0,TT), ylim = c(min(y), max(y)+max(y)/5),
     las = 1)
points(1:TT, data$GPP,
       pch = 19, cex = 0.7, col="red", ty = "o")
legend("top",
       legend = c("Obs.", "True states"),
       pch = c(3, 19),
       col = c("blue", "red"),
       lty = c(3, 1),
       horiz=TRUE, bty="n", cex=0.9)


TT<-1098
time=1:TT
plot.ts(data$GPP)

##Plot ACF and PACF
acf(y)
pacf(y)

##Difference
y_diff=diff(y,differences=1)

##Plot ACF and PACF of difference
acf(y_diff)
pacf(y_diff)

###Create a y with some missing data
y_miss<-GPP
y_na<- which(y %in% sample(y,20)) ##2,5,15
y_miss[y_na]<-NA
y_miss=as.data.frame(y_miss)

##Create objects with the location index of missing or observed data AND number
## of missing or observed data

missing <- lapply(y_miss,function(var){which(is.na(var))}) # Identify the rows for each variable in airquality with NA
observed <- lapply(y_miss,function(var){which(!is.na(var))}) # Identify the rows for each variable in airquality with observed data
nMissing <- lapply(missing,function(var){length(var)}) # How many NAs?
nObserved <- lapply(observed,function(var){length(var)}) # How many NAs?
names(missing) <- paste0(names(missing),'_index_mis') # Rename indices to variable_index
names(nMissing) <- paste0(names(nMissing),'_nMiss') # Rename n_miss to variable_nMiss
names(observed) <- paste0(names(observed),'_index_obs') # Rename indices to variable_index
names(nObserved) <- paste0(names(nObserved),'_nObs') # Rename n_Obs to variable_nObs


##Replace NAs with arbitrary number to make Stan happy
y_miss[is.na(y_miss)] <- -100


##Stan Code: State Space

sink("ss_reg_missing data.stan")

cat("
    data {
    int TT; // Latent state variable
    int N; // Number of observations
    int Y; // number of years
    vector[N] y_miss; // Response variable, including missing values
    int y_nMiss;
    int y_index_mis[y_nMiss];
   
    }
    
    parameters {
    vector[y_nMiss] y_imp;// Missing data
    vector[TT] X; // latent state data
    real<lower = 0> sdo; //  standard deviation observation
    real<lower = 0> sdp; //  standard deviation process
    real phi;  // auto-regressive parameter
   
    }
    
    transformed parameters { 
    vector[N] y; // makes the data a transformed variable
    y=y_miss; // 
    y[y_index_mis] =y_imp; // replaces missing data in y with estimated data
    } 
    
    model {
    for (j in 1:Y){
    
    X[1,j]<-y[1,j]  // initialization 
    mu[1,j]<-X[1,j] // initialization
    
    for (t in 2:TT){
    X[t,j] ~ normal(mu[t,j], sdp); // process model
    mu[t,j]<-phi[j]*X[i-1,j]+b0[j]+b1[j]*nit[i,j]+b2[j]*light[i,j]+b3[j]*srp[i,j] //regression model with AR errors
      
    y[t,j] ~ normal(X[t,j], sdo); // observation model
      
    // start of priors for each year
    sdp ~ normal(0, 2); 
    sdo ~ normal(0, 2);
    
    b0[j]~dnorm(b0mean,pow(b0sd,-2))
    b1[j]~dnorm(b1mean,pow(b1sd,-2))
    b2[j]~dnorm(b2mean,pow(b2sd,-2))
    b3[j]~dnorm(b3mean,pow(b3sd,-2))
    
    phi[j]~dbeta(a,b)T(0.001,0.999)
    
    // start of single parameters priors and hyperpriors
    obserr <- pow(sigmaobs, -2)
    sigmaobs ~ dunif(0, 1)
    b0mean~dnorm(0,pow(5,-2))
    b0sd~ dunif(0, 5)
    b1mean~dnorm(0,pow(5,-2))
    b1sd~dunif(0, 10)
    b2mean~dnorm(0,pow(5,-2))
    b2sd~dunif(0, 10)
    b3mean~dnorm(0,pow(5,-2))
    b3sd~dunif(0, 10)
    
    //Prior for AR coefficient needs to be reparameterized
    a<-meantheta*kappa 
    b<-(1-meantheta)*kappa
    meantheta~dbeta(1,1)
    kappa~dgamma(1,0.1)
    procerr <- pow(sigmaproc, -2)
    sigmaproc ~ dunif(0, 5)
    
    }
    }
    
    generated quantities{
    y.rep[i,j]~dnorm(mu[i,j], sdp)
    
    }
    "
    
    ,fill=TRUE)
sink()
closeAllConnections()

##Load data

data <- list(N = length(y_miss$y_miss),
                TT = length(y_miss$y_miss),
                y_nMiss = nMissing$y_miss_nMiss,
                y_nObs = nObserved$y_miss_nObs,
                y_index_mis = missing$y_miss_index_mis,
                y_index_obs = observed$y_miss_index_obs,
                y_miss= y_miss$y_miss
)

##Run Stan

fit <- stan("ts_ss_missing data.stan", data = data,  iter = 5000, chains = 4)

print(fit)
class(fit)
traceplot(fit, pars=c("sdo", "sdp","phi" ))
rstan::check_hmc_diagnostics(fit)
fit_extract<-rstan::extract(fit)

##Create object with estimated missing data

fit_summary<-summary(fit, probs=c(0.05,.5,.95))$summary
my_data <- as_tibble(fit_summary)
y_est_data<-my_data %>% slice(1:length(missing$y_miss_index_mis))

##Create object with observed data

y_miss_index<-missing$y_miss_index_mis
y_obs_data<-as.data.frame(cbind(time[missing$y_miss_index_mis],y[missing$y_miss_index_mis]))
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
raw.data=as.data.frame(cbind(time2, y_diff))

ggplot(data=y_combined, aes(x=date, y=estimated))+
  geom_line(data=raw.data, aes(x=time2, y=y_diff, color="black"), size=1)+
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

