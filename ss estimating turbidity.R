


##Turbidity--Duwamish River, Tukwila, Washington (https://www.sciencebase.gov/catalog/item/5a1dba7de4b09fc93dd7c022)
fake.turb<-read.csv("fake_turb.csv")

summary(fake.turb)

##Force some missing data
#Note-the model will still work if no data is missing but the objects in the next block of code (i.e., "y_index_mis", etc.)
#still need to be created
y_miss<-fake.turb$turb
#y_na<- which(y_miss %in% sample(y_miss,50)) ##the X in (y_miss,X) is the number of created missing data points 
#turb_miss[y_na]<-NA

##Create objects with the location (index) of missing or observed data AND number
## of missing or observed data
y_index_mis <- which(is.na(y_miss)) # Identify the rows for each variable in airquality with NA
y_index_obs <- which(!is.na(y_miss)) # Identify the rows for each variable in airquality with observed data
y_nMiss <- length(y_index_mis)# How many NAs?
y_nObs <- length(y_index_obs) # How many NAs?


##Replace NAs with arbitrary number to make Stan happy
y_miss[is.na(y_miss)] <- 10000



sink("ss estimate turb.stan")

cat("
    data {
    int<lower=1> N; // Number of observations
    vector[N] y_miss; // Response variable, including missing values
    int y_nMiss; // number of missing values
    int y_index_mis[y_nMiss]; // index or location of missing values within the dataset
    real Q[N]; //discharge data
    }
    
    parameters {
    vector<lower = 0>[y_nMiss] y_imp;// Missing data
    vector<lower = 0>[N] X; // latent state data
    real<lower = 0> sigma_proc; // process error
    real<lower = 0> sigma_obs; // observation error
    real<lower = 0, upper=1 > phi;  // auto-regressive parameter
    real b1; // discharge parameter
    }
    
    transformed parameters { 
    vector[N] y;
    y=y_miss; // makes the data a transformed variable
    y[y_index_mis] =y_imp; // replaces missing data in y with estimated data
    } 
    
    model {
    X[1]~lognormal(y[1],sigma_proc); //set initial state
    

    y ~ lognormal(X, sigma_obs); // observation model
    
    
    for (t in 2:N){
    X[t] ~ lognormal(phi*X[t-1]+b1*Q[t], sigma_proc); // process model with unknown process error //regression model with AR errors
    
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
 x_rep[1]=lognormal_rng(y[1], 0.1);
 
 for (t in 2:N) {
 x_rep[t]=lognormal_rng(phi*X[t-1]+b1*Q[t], sigma_proc);
 }
 for (t in 1:N) {
 y_rep[t]=lognormal_rng(x_rep[t], sigma_obs);
 }

  
 }
 
    "
    ,fill=TRUE)
sink()
closeAllConnections()


data <- list(   N = length(y_miss),
                y_nMiss = y_nMiss,
                y_nObs = y_nObs,
                y_index_mis =y_index_mis,
                y_index_obs = y_index_obs,
                y_miss= y_miss,
                Q=fake.turb$Q
                
)

fit.turb<- stan("ss estimate turb.stan", data = data,  iter = 2000, chains = 3)

##Print and extract
fit_extract<-rstan::extract(fit.turb)
print(fit.turb, pars=c( "sigma_proc","phi", "b1", "sigma_obs" ))

##HMC diagnostics
rstan::check_hmc_diagnostics(fit.turb)

##Traceplots
traceplot(fit.turb, pars=c( "sigma_proc", "sigma_obs" ))
traceplot(fit.turb, pars=c("phi", "b1"))

##Posterior densities comopared to known parameters
plot_sdo <- stan_dens(fit.turb, pars="sigma_obs") 
plot_sdp <- stan_dens(fit.turb, pars="sigma_proc") 
plot_phi <- stan_dens(fit.turb, pars="phi") 
plot_b1 <- stan_dens(fit.turb, pars="b1")+ xlab("Q beta")
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

fit_summary<-summary(fit.turb, probs=c(0.05,.5,.95))$summary

my_data <- as_tibble(fit_summary)
y_est_data<-my_data %>% slice(1:length(turb_index_mis))

##Create object with observed data

turb_miss_index<-data$turb_index_mis
date<-fake.turb$date
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





