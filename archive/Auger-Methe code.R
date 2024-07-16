# Create a vector that will keep track of the states
# It's of length T + 1 (+1 for t=0)
# T is not a good name in R, because of T/F, so we use TT
TT <- 700
t<-1:TT
z <- numeric(TT)
# Standard deviation of the process variation
sdp <- 0.1
phi <-0.8
b0<-0.1
b1<-0.1
z[1]<-3.5
# Set the seed, so we can reproduce the results
set.seed(553)
# For-loop that simulates the state through time, using i instead of t,
for (t in 2:TT){
  z[t] = b0+phi*z[t-1]+light.l[t]*b1+rnorm(1, 0, sdp)
}
  # Note that this index is shifted compared to equation in text,
  # because we assume the first value to be at time 0



plot(1:TT, z,
     pch = 19, cex = 0.7, col="red", ty = "o",
     xlab = "t", ylab = expression(z[t]), las=1)

# Create a vector that will keep track of the observations
# It's of length T
y <- numeric(TT)
# Standard deviation of the observation error
sdo <- 0.2
# For t=1, ... T, add measurement error
# Remember that z[1] is t=0
y <- z[2:(TT+1)] + rnorm(TT, 0, sdo)

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


##Force some missing data
#Note-the model will still work if no data is missing but the objects in the next block of code (i.e., "y_index_mis", etc.)
#still need to be created
set.seed(620)
y_miss<-z
y_na<- which(y_miss %in% sample(y_miss,350)) ##the X in (y_miss,X) is the number of created missing data points 
y_miss[y_na]<-NA
y_miss[1]<-z[1]

##Create objects with the location (index) of missing or observed data AND number
## of missing or observed data
y_index_mis <- which(is.na(y_miss)) # Identify the rows for each variable in airquality with NA
y_index_obs <- which(!is.na(y_miss)) # Identify the rows for each variable in airquality with observed data
y_nMiss <- length(y_index_mis)# How many NAs?
y_nObs <- length(y_index_obs) # How many NAs?

##Replace NAs with arbitrary number to make Stan happy
y_miss[is.na(y_miss)] <- -100



sink("toy2p.stan")

cat("
 
/*----------------------- Data --------------------------*/
  /* Data block: defines the objects that will be inputted as data */
  data {
    int TT; // Length of state and observation time series
    vector[TT] y_miss; // Observations
    real z0; // Initial state value
    vector[TT] light;
    int y_nMiss; // number of missing values
    int y_index_mis[y_nMiss]; // index or location of missing values within the dataset
  }
/*----------------------- Parameters --------------------------*/
  /* Parameter block: defines the variables that will be sampled */
  parameters {
    vector[y_nMiss] y_imp;// Missing data
    real<lower=0> sdp; // Standard deviation of the process equation
    //real<lower=0> sdo; // Standard deviation of the observation equation
    real b0;
    real b1;
    real<lower = 0, upper=1 > phi; // Auto-regressive parameter
    //vector[TT] z; // State time series
  }
  transformed parameters { 
    vector[TT] y;
    y=y_miss; // makes the data a transformed variable
    y[y_index_mis] =y_imp; // replaces missing data in y with estimated data
    } 
  /*----------------------- Model --------------------------*/
  /* Model block: defines the model */
  model {
    // Prior distributions
    //sdo ??? normal(0, 1);
    sdp ??? cauchy(0, 1);
    phi ~ beta(1,1);
    b0 ~ normal(0,5);
    b1 ~ normal(0,5);
    // Distribution for the first state
    y[1] ??? normal(z0, sdp);
    // Distributions for all other states
    for(t in 2:TT){
      y[t] ??? normal(b0+y[t-1]*phi+light[t]*b1, sdp);
    }
    // Distributions for the observations
   //for(t in 1:TT){
    
    //  y[t] ??? normal(z[t], sdo);
   //}
  }
    "
    ,fill=TRUE)
sink()
closeAllConnections()

dataStan <- list(y_miss=y_miss[[3]], 
                 TT=N, 
                 z0=z[1], 
                 light=light.l,
                 y_nMiss = y_nMiss[[3]],
                 y_index_mis =y_index_mis[[3]])

f2pStan <- stan(file = "toy2p.stan",
                data = dataStan,
                chains = 4, iter = 4000)

traceplot(f2pStan, pars=c("sdp", "phi", "b0" ,"b1"))

print(f2pStan, max = 50)

##Posterior densities comopared to known parameters
#plot_sdo <- stan_dens(f2pStan, pars="sigma_obs") + geom_vline(xintercept =sd.o)
plot_sdp <- stan_dens(f2pStan, pars="sdp") + geom_vline(xintercept = sdp)+xlab("Error")+geom_vline(xintercept = 0.104, col="blue")+geom_vline(xintercept = 0.10250743, col="yellow")
plot_phi <- stan_dens(f2pStan, pars="phi") + geom_vline(xintercept = phi)+xlab("AR beta")+ geom_vline(xintercept = 0.802, col="blue")+ geom_vline(xintercept = 0.81164313, col="yellow")
plot_b1 <- stan_dens(f2pStan, pars="b1") + geom_vline(xintercept = b1)+xlab("Light beta")+geom_vline(xintercept = 0.0917, col="blue")+geom_vline(xintercept = 0.08670065, col="yellow")
grid.arrange(plot_b1, plot_phi,plot_sdp, nrow=1)

##Create object with estimated missing data

fit_summary<-summary(f2pStan, probs=c(0.05,.5,.95), pars=c("sdp", "phi", "b0" ,"b1"))$summary

my_data <- as_tibble(fit_summary)
y_est_data<-my_data %>% slice(1:length(y_index_mis[[1]]))

##Create object with observed data

y_miss_index<-dataStan$y_index_mis
date<-1:N
y_obs_data<-as.data.frame(cbind(date[dataStan$y_index_mis],x[dataStan$y_index_mis]))
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
  geom_errorbar(aes(ymin=low, ymax=high), width=0.01)+
  theme(legend.position="top")+
  ylab("Estinated data (+/- 95% Cred. Interval)")+
  xlab("Simulated data")+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))

##Plot observed and estimated time series
t<-1:N
raw.data=as.data.frame(cbind(t, x))

ggplot(data=y_combined, aes(x=date, y=estimated))+
  geom_point(data=raw.data, aes(x=t, y=x, color="gray"), size=3)+
  geom_point(aes(color="red"), size=3)+
  geom_errorbar(data=y_combined,aes(x=date,ymin=low, ymax=high))+ 
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  scale_color_identity(guide = "legend", name=element_blank(), labels=c("Raw data", "Estimated points"))+
  theme(legend.position="top")+
  ylab("Data")+
  xlab("Time")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


