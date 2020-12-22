
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
dat <- rnorm_multi(n = 365,
                   mu = c(0,0,0,0),
                   sd = c(5,5,0.5,0.5),
                   r = c(1,.6,.1,.3,.6,1,.5,.5,.1,.5,1,.9,.3,.5,.9,1), ##correlation matrix
                   varnames = c("srp", "nit", "light"),
                   empirical = FALSE)
dat1$Y<-rep(1,366)
dat1$t<-seq(from=1, to=366, by=1)

TT=365

Y<-numeric(TT+1)

b0<-1
b1<-.8
sd.p<-.2
for (t in 1:TT){
Y[t+1] = .9*Y[t]+b0+b1*dat1$light[t]+rnorm(1, 0, sd.p)
}

plot(0:TT, Y,
     pch = 19, cex = 0.7, col="red", ty = "o",
     xlab = "dat1$t", ylab = expression(z[t]), las=1)
plot(Y~dat1$nit)

quad<-abs(seq(-2,2, length=366))
dat1$light<-dat1$light-quad
dat1$GPP<-dat1$GPP-quad
plot(dat1$t, dat1$GPP)
plot(dat1$GPP~dat1$light)

dat2 <- rnorm_multi(n = 366,
                    mu = c(0,0,0,0),
                    sd = c(5,5,0.5,0.5),
                    r = c(1,.4,.1,.3,.4,1,.5,.8,.1,.5,1,.9,.3,.8,.9,1), ##correlation matrix
                    varnames = c("srp", "nit", "light", "GPP"),
                    empirical = FALSE)
dat2$Y<-rep(2,366)
dat2$t<-seq(from=1, to=366, by=1)
quad<-abs(seq(-2,2, length=366))
dat2$light<-dat2$light-quad
dat2$GPP<-dat2$GPP-quad
plot(dat2$t, dat2$light)


dat3 <- rnorm_multi(n = 366,
                    mu = c(0,0,0,0),
                    sd = c(5,5,0.5,0.5),
                    r = c(1,.4,.1,.3,.4,1,.5,.8,.1,.5,1,.6,.3,.8,.6,1), ##correlation matrix
                    varnames = c("srp", "nit", "light"),
                    empirical = FALSE)

dat3$Y<-rep(3,366)
dat3$t<-seq(from=1, to=366, by=1)
quad<-abs(seq(-5,5, length=366))
dat3$light<-dat3$light-quad
dat3$GPP<-dat3$GPP-quad


data<-rbind(dat1,dat2,dat3)





###Create observed data
# Create a vector that will keep track of the observations
# It's of length T
TT=1098
y <- numeric(TT)
# Standard deviation of the observation error
sdo <- 1
# For t=1, ... T, add measurement error
# Remember that z[1] is t=0
y <- Y[2:(TT+1)] + rnorm(TT, 0, sdo)



###Plot simulated and observed data
plot(1:TT, Y[2:366],
     pch=3, cex = 0.8, col="blue", ty="o", lty = 3,
     xlab = "t", ylab = expression(z[t]),
     #xlim = c(0,TT), ylim = c(min(y), max(y)+max(y)/5),
     las = 1)
points(1:TT, y,
       pch = 19, cex = 0.7, col="red", ty = "o")
legend("top",
       legend = c("Obs.", "True states"),
       pch = c(3, 19),
       col = c("blue", "red"),
       lty = c(3, 1),
       horiz=TRUE, bty="n", cex=0.9)

sd(y)

summary (y)
##Plot ACF and PACF
acf(y)
pacf(y)

##Difference
y_diff=diff(y,differences=1)

##Plot ACF and PACF of difference
acf(y_diff)
pacf(y_diff)

###Create a y with some missing data
y_miss<-y
y_na<- which(y_miss %in% sample(y_miss,100)) ##2,5,15
y_miss[y_na]<-NA
#y_miss=as.data.frame(y_miss)

##Create objects with the location index of missing or observed data AND number
## of missing or observed data

y_index_mis <- which(is.na(y_miss)) # Identify the rows for each variable in airquality with NA
y_index_obs <- which(!is.na(y_miss)) # Identify the rows for each variable in airquality with observed data
y_nMiss <- length(y_index_mis)# How many NAs?
y_nObs <- length(y_index_obs) # How many NAs?



##Replace NAs with arbitrary number to make Stan happy
y_miss[is.na(y_miss)] <- -100


##Stan Code: State Space

sink("ss_reg_missing data.stan")

cat("
    data {
    int TT; // Number of observations
    vector[TT] y_miss; // Response variable, including missing values
    int y_nMiss; // number of missing values
    int y_index_mis[y_nMiss]; // index or location of missing values within the dataset
    real nit[TT];
    real light[TT];
    real srp[TT];
    }
    
    parameters {
    vector[y_nMiss] y_imp;// Missing data
    vector[TT] X; // latent state data
    real<lower = 0> sdo; //  standard deviation observation
    real<lower = 0> sdp; //  standard deviation process
    real phi;  // auto-regressive parameter
    real b0; // intercept
    real b1; // nit parameter 
    real b2; // light parameter
    real b3; // srp parameter 
    }
    
    transformed parameters { 
    vector[TT] y;
    y=y_miss; // makes the data a transformed variable
    y[y_index_mis] =y_imp; // replaces missing data in y with estimated data
    } 
    
    model {
    
    X[1]~normal(0,2);
    y[1]~normal(0,2);
    
    for (t in 2:TT){
    X[t] ~ normal(phi*X[t-1]+b0+b1*nit[t]+b2*light[t]+b3*srp[t], sdp); // process model //regression model with AR errors
    
    y[t] ~ normal(X[t], sdo); // observation model
    }  
    // start of priors for each year
    sdp ~ normal(0, 2); 
    sdo ~ normal(0, 2);
    
    // start of single parameters priors and hyperpriors
    
    b0~normal(0, 5);
    b1~normal(0, 10);
    b2~normal(0, 10);
    b3~normal(0, 10);
    
    //Prior for AR coefficient needs to be reparameterized
    phi~normal(0,5);
    }
    
    "
    
    ,fill=TRUE)
sink()
closeAllConnections()

##Load data

data <- list(   TT = length(y_miss),
                y_nMiss = y_nMiss,
                y_nObs = y_nObs,
                y_index_mis =y_index_mis,
                y_index_obs = y_index_obs,
                y_miss= y_miss,
                nit=dat1$nit,
                srp=dat1$srp,
                light=dat1$light
)
summary(data)
##Run Stan

fit <- stan("ss_reg_missing data.stan", data = data,  iter = 5000, chains = 3)

print(fit)
class(fit)
traceplot(fit, pars=c("sdo", "sdp","phi" ))
traceplot(fit, pars=c("b0", "b1","b2","b3" ))
rstan::check_hmc_diagnostics(fit)
fit_extract<-rstan::extract(fit)

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

