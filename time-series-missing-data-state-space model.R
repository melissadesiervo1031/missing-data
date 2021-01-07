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
set.seed(620)

###Simulate data
TT=365 #length of data
t=0:TT #time

##light
light <- rnorm(n = TT, mean = 0, sd = 0.5)

##GPP based on time-series model with known paramters
x=length(TT) #for storage
b0<-1 #intercept
phi<--0.8
b1<-1.5 #light
sd.p<-1 #process error

#Estimate latent state
for (t in 1:TT){
        x[t+1] = phi*x[t]+b0+b1*light[t]+rnorm(1, 0, sd.p)
}
#Estimate observed data
sdo <- rnorm(TT, 0, 1.5)
y <- x[2:(TT+1)] + sdo

#Create fixed observation for inclusion in  model
sdo_sd<-abs((sdo))


##Plot observed and latent
plot(1:TT, x[2:366],
     pch=3, cex = 0.8, col="blue", ty="o", lty = 3,
     xlab = "t", ylab = expression(z[t]),
     xlim = c(0,TT), ylim = c(min(y), max(y)+max(y)/5),
     las = 1)
points(1:TT, y,
       pch = 19, cex = 0.7, col="red", ty = "o")
legend("top",
       legend = c("Obs.", "Latent states"),
       pch = c(3, 19),
       col = c("blue", "red"),
       lty = c(3, 1),
       horiz=TRUE, bty="n", cex=0.9)

##Force some missing data
#Note-the model will still work if no data is missing but the objects in the next block of code (i.e., "y_index_mis", etc.)
#still need to be created
y_miss<-y
#y_na<- which(y_miss %in% sample(y_miss,50)) ##the X in (y_miss,X) is the number of created missing data points 
#y_miss[y_na]<-NA

##Create objects with the location (index) of missing or observed data AND number
## of missing or observed data
y_index_mis <- which(is.na(y_miss)) # Identify the rows for each variable in airquality with NA
y_index_obs <- which(!is.na(y_miss)) # Identify the rows for each variable in airquality with observed data
y_nMiss <- length(y_index_mis)# How many NAs?
y_nObs <- length(y_index_obs) # How many NAs?

##Replace NAs with arbitrary number to make Stan happy
y_miss[is.na(y_miss)] <- -100

##Stan Code

sink("ss_reg_missing data.stan")

cat("
    data {
    int TT; // Number of observations
    vector[TT] y_miss; // Response variable, including missing values
    int y_nMiss; // number of missing values
    int y_index_mis[y_nMiss]; // index or location of missing values within the dataset
    real light[TT]; //light data
    real sdo[TT]; //fixed observation error sd
    }
    
    parameters {
    vector[y_nMiss] y_imp;// Missing data
    vector[TT] X; // latent state data
    real<lower = 0> sdp; // process error
    real phi;  // auto-regressive parameter
    real b0; // intercept
    real b1; // light parameter 
    }
    
    transformed parameters { 
    vector[TT] y;
    y=y_miss; // makes the data a transformed variable
    y[y_index_mis] =y_imp; // replaces missing data in y with estimated data
    } 
    
    model {
    for (t in 2:TT){
    X[t] ~ normal(phi*X[t-1]+b0+b1*light[t], sdp); // process model with unknown process error //regression model with AR errors
    
    y[t] ~ normal(X[t], sdo[t]); // observation model with fixed observation error
        }  
    
    // error priors
    sdp ~ cauchy(0,1); 
  
    
    // single parameters priors 
    b0~normal(0, 1);
    b1~normal(0, 1);
   
    //Prior for AR coefficient needs to be reparameterized
    phi~normal(0,1);
    }
    
    generated quantities {
    vector[TT] log_y_rep; // replications from posterior predictive dist

    for (t in 2:TT) {
     for (n in 1:TT) {
    real log_y_hat_n=phi*X[t-1]+b0+b1*light[t];
    log_y_rep[n]=normal_rng(log_y_hat_n, sdp);
    }
    }
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
                light=light,
                sdo=sdo_sd
)


##Run Stan

fit.sdo <- stan("ss_reg_missing data.stan", data = data,  iter = 5000, chains = 3)

print(fit)
rstan::check_hmc_diagnostics(fit.sdo)
fit_extract<-rstan::extract(fit.sdo)
traceplot(fit.sdo, pars=c( "sdp","phi" ))
traceplot(fit.sdo, pars=c("b0", "b1"))

#pairs(fit, pars = c("sdo", "sdp", "b0","phi","b1","lp__"), las = 1)

plot_sdo <- stan_dens(fit.sdo, pars="sdo") + geom_vline(xintercept =sdo)
plot_sdp <- stan_dens(fit.sdo, pars="sdp") + geom_vline(xintercept = sd.p)
plot_phi <- stan_dens(fit.sdo, pars="phi") + geom_vline(xintercept = phi)
plot_b1 <- stan_dens(fit.sdo, pars="b1") + geom_vline(xintercept = b1)
plot_b0 <- stan_dens(fit.sdo, pars="b0") + geom_vline(xintercept = b0)
grid.arrange(plot_b0,plot_b1, plot_phi,plot_sdp, nrow=2)



##Posterioir Predictive check
yrep1 <- fit_extract$log_y_rep
samp100 <- sample(nrow(yrep1), 100)
ppc_dens_overlay(y_miss, yrep1[samp100, ]) 
ppc_stat(y_miss, yrep1[samp100, ])

##Manual calculation of posterior predictive check (Did I do it right in Stan?)
# extract the samples
nsims<-dim(fit_extract$sdp)[1]

X_s<-fit_extract$X
sdp_s<-fit_extract$sdp
phi_s<-fit_extract$phi
b0_s<-fit_extract$b0
b1_s<-fit_extract$b1

# empty matrix to store samples
y_rep <- matrix(NA, nrow = TT, ncol = nsims)

for (t in 2:TT){
  y_rep[t,] = rnorm(nsims, mean=phi_s*X_s[t-1]+b0_s+b1_s*light[t],sd= sdp_s)
}  

colnames(y_rep) <- 1:nsims

dr <- as_tibble(y_rep)
dr <- dr %>% bind_cols(i = 1:TT, y_obs = y_miss)

dr <- dr %>% 
  pivot_longer(`1`:`100`, names_to = "sim", values_to = "y_rep")

data.df<- data.frame(matrix(unlist(data), nrow=365, byrow=T),stringsAsFactors=FALSE)

dr %>% 
  filter(dr %in% samp100) %>% 
  ggplot(aes(y_rep, group = sim)) + 
  geom_density(alpha = 0.2, aes(color = "y_rep")) + 
  geom_density( data=data.df%>% mutate(sim = 1), 
               aes(x = y_miss, col = "y")) + 
  scale_color_manual(name = "", 
                     values = c("y" = "darkblue", 
                                "y_rep" = "lightblue")) + 
  theme_bw(base_size = 16)




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





