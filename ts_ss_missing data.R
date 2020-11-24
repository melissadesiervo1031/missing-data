
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(dplyr)
library(shinystan)

###Simulate data
# Create a vector that will keep track of the states
# It's of length T + 1 (+1 for t=0)
# T is not a good name in R, because of T/F, so we use TT
TT <- 200
time<-0:TT
z <- numeric(TT + 1)
# Standard deviation of the process variation
sdp <- 0.5
# Set the seed, so we can reproduce the results
set.seed(620)
# For-loop that simulates the state through time, using i instead of t,
for(i in 1:TT){
  # This is the process equation
  z[i+1] <- z[i] + rnorm(1, 0, sdp)
 }

###Plot simulated data
plot(0:TT, z,
     pch = 19, cex = 0.7, col="red", ty = "o",
     xlab = "t", ylab = expression(z[t]), las=1)


###Create observed data
# Create a vector that will keep track of the observations
# It's of length T
y <- numeric(TT)
# Standard deviation of the observation error
sdo <- .5
# For t=1, ... T, add measurement error
# Remember that z[1] is t=0
y <- z[2:(TT+1)] + rnorm(TT, 0, sdo)


###Plot simulated and observed data
plot(1:TT, y,
     pch=3, cex = 0.8, col="blue", ty="o", lty = 3,
     xlab = "t", ylab = expression(z[t]),
     xlim = c(0,TT), ylim = c(min(y), max(y)+max(y)/5),
     las = 1)
points(0:TT, z,
       pch = 19, cex = 0.7, col="red", ty = "o")
legend("top",
       legend = c("Obs.", "True states"),
       pch = c(3, 19),
       col = c("blue", "red"),
       lty = c(3, 1),
       horiz=TRUE, bty="n", cex=0.9)

##Plot ACF and PACF
acf(y)
pacf(y)

##Difference
y_diff=diff(y,differences=1)

##Plot ACF and PACF of difference
acf(y_diff)
pacf(y_diff)

###Create a y with some missing data
y_miss<-y_diff
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

sink("ts_ss_missing data.stan")

cat("
    data {
    int TT; //Latent state variable
    int N; //Number of observations
    vector[N] y_miss; //Response variable, including missing values
    int y_nMiss;
    int y_index_mis[y_nMiss];
   
    }
    
    parameters {
    vector[y_nMiss] y_imp;//Missing data
    vector[TT] x; // latent state data
    real<lower = 0> sdo; //  standard deviation observation
    real<lower = 0> sdp; //  standard deviation process
    real phi;  // auto-regressive coefficient
   
    }
    
    transformed parameters { 
    vector[N] y;
    y=y_miss;
    y[y_index_mis] =y_imp;
    } 
    
    model {
    for (t in 2:TT){
    x[t] ~ normal(phi*x[t-1], sdp);
      }
     for (t in 1:TT){
    y[t] ~ normal(x[t], sdo);
      }
    sdp ~ normal(0, 2);
    sdo ~ normal(0, 2);
   
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

fit <- stan("ts_ss_missing data.stan", data = data,  iter = 3000, chains = 4)

print(fit)
class(fit)
traceplot(fit, pars=c("sdo", "sdp","phi" ))
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

