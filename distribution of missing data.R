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
library(tidyr)
library(MASS)
library(plyr)
library(ggplot2)
library(countreg)
library(abind)
# SDW time series ---------------------------------------------------------


#Load Data
ts.SDW<-read.csv("SDW.ts_with light.csv")
ts.SDW$Year<-as.factor(ts.SDW$year)

##Format date for filling in dates with missing data
ts.SDW$date<-as.Date(ts.SDW$date,format="%m/%d/%Y")

##Remove negative GPP values
ts.SDW<- subset(ts.SDW,GPP >= 0)

##Fill in missing dates
df<-complete(date = seq.Date(min(ts.SDW$date), max(ts.SDW$date), by="day"), data=ts.SDW)

length(which(is.na(df$GPP), TRUE))/length(df$GPP)

##Code data to be missing=1 or observed=0
N<-ifelse(is.na(df$GPP),1,0)

## Calculate length of data gaps
#prepare vectors for loop
x<-NA 
length_miss<-NA
x[1]<-N[1]
length_miss[1]<-1

#loop over data and return the length of missing data gaps
for (i in 2:length(df$GPP)){
  x[i]<-ifelse(N[i]==1,1+x[i-1], 0) #If data is missing add to yesterdays value by 1. If data observed then 0.
  length_miss[i]<-ifelse(x[i]==0,x[i-1], NA) #Save length of gap on the day data is observed (0) after a gap
 # length_miss[i]<-ifelse(length_miss[i]==0,NA,length_miss[i]) #Replace days with observed data with NA
}

## Calculate length of observed data before gap. Same as above but for lengths of observed data
z<-NA
z[1]<-N[1]
length_obs<-NA
for (i in 2:length(df$GPP)){
  z[i]<-ifelse(N[i]==0,1+z[i-1], 0)
  length_obs[i]<-ifelse(z[i]==0,z[i-1], NA)
 # length_obs[i]<-ifelse(length_obs[i]==0,NA,length_obs[i] )
}


## save to full dataset and check that everything is correct
df$length_obs<-length_obs
df$length_miss<-length_miss

## estimate negative binomial distribution parameters of missing data length integers
length_miss_parm<-length_miss[complete.cases(length_miss)]
length_miss_parm_zero<-length_miss_parm[length_miss_parm<90]
length_miss_parm<-length_miss_parm_zero[length_miss_parm_zero>0]

##regular neg binomial
func<- function(P){
  -sum( dnbinom (length_miss_parm, size= P[1], prob = P[2], p0=P[3]))
}

mle<-nlm(func,P<- c(2.5,0.5,0.5))


##zero truncated negative binomal
func<- function(P){
  -sum( dzmnbinom (length_miss_parm_zero, size= P[1], prob = P[2], p0=0.9), log=TRUE)
}

mle<-nlm(func,P<- c(2,0.5))



## plot histogram
hist(length_miss_parm_zero, breaks=seq(0,max(length_miss[complete.cases(length_miss)]),by=1), prob=TRUE, xlim=c(0,70),ylim=c(0,1), xlab="Length of missing data gap",
     #main="Negative binomial distribution with mu=8.6 and size=0.65")
     main="Histogram of missing data gap lengths from Shatto")

## plot distribution of negative binomial using estimated parameters
x=0:max(length_miss_parm)
dist<-dzmnbinom(x=x,size=mle$estimate[1],prob=mle$estimate[2], p0=0.9)
points(dist, col="red")
plot(x, dist)

## create random missing data vector using zero-inflated distribution
sim.prob<-rzmnbinom(x, size=mle$estimate[1],size=25,p0=0.8)
hist(sim.prob, prob=TRUE,ylim=c(0,1) )


## create random missing data vector using by randomly pulling from the data
N<-1:365
n<-365
s<-sample(length_miss_parm_zero,300,replace=TRUE)
index<-NA

obs<-na.omit(length_obs)
obs<-obs[which(obs>0)]
miss<-na.omit(length_miss)
miss<-miss[which(miss>0)]

saveRDS(miss, file = "SDW_miss_length.RDS")
saveRDS(obs, file = "SDW_obs_length.RDS")


sample.miss<-sample(miss,200,replace=TRUE)
sample.obs<-sample(obs,200,replace=TRUE)

N<-1:365
n<-365
new<-sample.obs[1]


for(i in 1:n){

remove<-seq(from=new+1, to=new+sample.miss[i], by=1)
N[remove]<-0
new<-new+sample.miss[i]+sample.obs[i+1]
if (new>=n){
  break
}
}


for (i in 1:N){
  s[i]<-sample(length_miss_parm_zero,1,replace=TRUE)

if (s[i]==0){
  new_N[i]<-N[i]
}else{
  new_N[N[i]:(N[i]+s[i])]<-0
}
  if (length(new_n)>= N)  {
    break
}

}


## plot sum of gaps by julian day (the day the gap ended)
#No observable trend. Suggests there is no bias of missing data for a specific time of the year
#Missing data seems to be randomly dispersed throughout the year
df %>% 
  group_by(jd) %>% 
  summarise(total = sum(length_miss,na.rm=TRUE)) %>% 
  with(plot(jd, total, ylab="Total missing data by Julian Day", xlab="Julian day"))

# Kalamath River (SV) time series ---------------------------------------------------------
setwd("C:/Users/matt/Documents/Laurel/SV_metab")

#Load Data
ts.SV<-read.csv("params.csv")
ts.SV$year<-(year(ts.SV$date))

##Format date for filling in dates with missing data
ts.SV$date<-as.Date(ts.SV$date,format="%Y-%m-%d")

##Remove negative GPP values
ts.SV<- subset(ts.SV,GPP.daily >= 0)
ts.SV<- subset(ts.SV, year>=2018)
ts.SV$year<-as.factor(year(ts.SV$date))

##Fill in missing dates
df<-complete(date = seq.Date(min(ts.SV$date), max(ts.SV$date), by="day"), data=ts.SV)

length(which(is.na(df$GPP.daily), TRUE))/length(df$GPP.daily)

##Code data to be missing=1 or observed=0
N<-ifelse(is.na(df$GPP.daily),1,0)

## Calculate length of data gaps
#prepare vectors for loop
x<-NA
length_miss<-NA
x[1]<-N[1]
length_miss[1]<-1

#loop over data and return the length of gaps
for (i in 2:length(df$GPP.daily)){
  x[i]<-ifelse(N[i]==1,1+x[i-1], 0) #If data is missing add to yesterdays value by 1. If data observed then 0.
  length_miss[i]<-ifelse(x[i]==0,x[i-1], NA) #Save length of gap on the day data is observed (0) after a gap
  length_miss[i]<-ifelse(length_miss[i]==0,NA,length_miss[i]) #Replace days with observed data with NA
}

## Calculate length of observed data before gap. Same as above but for lengths of observed data
z<-NA
z[1]<-N[1]
length_obs<-NA
for (i in 2:length(df$GPP.daily)){
  z[i]<-ifelse(N[i]==0,1+z[i-1], 0)
  length_obs[i]<-ifelse(z[i]==0,z[i-1], NA)
  length_obs[i]<-ifelse(length_obs[i]==0,NA,length_obs[i] )
}


## save to full dataset and check that everything is correct
df$length_obs<-length_obs
df$length_miss<-length_miss

## estimate negative binomial distribution parameters of missing data length integers
length_miss_parm<-length_miss[complete.cases(length_miss)]
length_miss_attempted<-sum(length_miss_parm[length_miss_parm<90])/sum(length_miss_parm)


## plot histogram
hist(length_miss_parm, breaks=seq(0,max(length_miss_parm+1),by=1), prob=TRUE, ylim=c(0,0.4), xlab="Length of missing data gap",
     #main="Negative binomial distribution with mu=8.6 and size=0.65")
     main="Histogram of missing data gap lengths from Shatto")

## plot distribution of negative binomial using estimated parameters
x=0:max(length_miss_parm)
dist<-dnbinom(x=x,size=ff[[1]][1],mu=ff[[1]][2])
points(dist, col="red")

## plot sum of gaps by julian day (the day the gap ended)
#No observable trend. Suggests there is no bias of missing data for a specific time of the year
#Missing data seems to be randomly dispersed throughout the year
df %>% 
  group_by(jd) %>% 
  summarise(total = sum(length_miss,na.rm=TRUE)) %>% 
  with(plot(jd, total, ylab="Total missing data by Julian Day", xlab="Julian day"))

# Appling et al. data ------------------------------------------------------

###Load and prep data
##MTT desktop
#setwd("C:/Users/mtrentman/IDrive-Sync/Postdoc/Estimating missing data/daily_predictions")
##Hall laptop
setwd("C:/Users/matt/IDrive-Sync/Postdoc/Estimating missing data/daily_predictions")
#sp<-read.table('daily.predictions.filled.tsv',header = TRUE) ##Missing dates filled in
sp<-read.csv('daily.predictions.filled.csv',header = TRUE) ##Missing dates filled in
sp$date.f<-as.Date(sp$date.f,format="%Y-%m-%d")
sp$site_name<-as.character(sp$site_name)
#sp<-read.table(file = 'daily_predictions.tsv',header = TRUE) ##Raw data
sd<-read.table(file = 'site_data.tsv',header = TRUE)##Site data
sd$site_name<-as.character(sd$site_name)

# Filling missing dates ---------------------------------------------------

#Format date for filling in dates with missing data
sp$date.f<-as.Date(sp$date,format="%Y-%m-%d")
sp$site_name<-as.factor(sp$site_name)
sp$year<-as.factor(format(as.Date(sp$date.f, format="%Y-%m-%d"),"%Y"))#Add year

#Remove negative GPP values
sp<- subset(sp,GPP >= 0)

###Fill in missing dates by site-year
t<-c()
sites<-(levels(sp[, 1]))
dummy1<-c()
for (i in 1:length(sites)){
  # for (g in levels(sp[,31])){
  dummy1 <-subset(sp,site_name==sites[i]) 
  #dummy2 <- sp[sp$year==g,]
  t<-rbind(t,complete(date.f = seq.Date(min(dummy1$date.f), max(dummy1$date.f), by="day"),data=dummy1))
}
#}

#Fill in site_name gaps
t<-as.data.frame(t %>% fill(site_name))
t$year<-as.factor(format(as.Date(t$date.f, format="%Y-%m-%d"),"%Y"))#Add year

write.table(t, "daily.predictions.filled.tsv",quote=FALSE, sep='\t',row.names = FALSE)


##I filled in and saved the data


# Filter sites ------------------------------------------------------------


# Filter sites
###Count NAs
cdata <- ddply(sp, c("site_name", "year"), summarize,
               N    = length(GPP),
               N.miss= sum(is.na(GPP)),
               Prop.miss=round(sum(is.na(GPP))/length(GPP)*100))


##Find sites with low missing data
full.year<-subset(cdata, N==365)
full.year<-subset(full.year,Prop.miss<90 )
###Summary of NAs               
#Histogram
hist(full.year$Prop.miss, main="Distribution of percent missing data from Appling et al. site*years", xlab="Percent missing data", cex=2)

med.miss<-median(full.year$Prop.miss)

appling<-ggplot(full.year, aes(x=Prop.miss)) + 
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, 
  color="#999999", fill="#999999", binwidth = 5)+
    theme(legend.position="top")+
  geom_vline(data=full.year,aes(xintercept=median(full.year$Prop.miss)), color="#FF00FF", size=2, linetype="dashed")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  theme(legend.position="top")+
  ylab("Density")+
  xlab("Percent missing data")+
  scale_x_continuous(limits=c(-5,105), breaks=c(0,20,40,60,80,100))+
  ggtitle("Distribution of percent missing data \n from Appling et al. 2018 with full site-years")+
  theme(plot.title = element_text(size = 15, face = "bold",hjust = 0.5))

ggsave("appling.tiff", width =8, height = 6)

high.miss.example=subset(sp, site_name=="nwis_01388000" & year==2010)
high_missing_integers<-which(is.na(high.miss.example$GPP))
high.miss.example$miss<-NA
high.miss.example$miss[high_missing_integers]<-5

ggplot(data=high.miss.example, aes(x=as.Date(date.f), y=GPP))+
  geom_point(color="black", size=3)+
  geom_point(aes(x=as.Date(date.f), y= miss),color="red", size=1)+
  geom_errorbar(aes(x=as.Date(date.f),ymin=GPP.lower, ymax=GPP.upper), width=0.2, size=0.5)+
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  theme(legend.position="top")+
  ylab("GPP")+
  xlab("Time (days)")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_x_date(limits = c(ymd(2013-01-01),ymd(2013-12-31)))+
  theme(plot.title = element_text(size = 15, face = "bold",hjust = 0.5))+
  ggtitle(paste0(sd$long_name[sd$site_name=="nwis_01388000"]," ","(2010)","
    ","25% missing"))

ggsave("low.miss.example.tiff", width =8, height = 6)

#summary stats
summary(full.year$Prop.miss)

##Find sites with low missing data
low.miss<-cdata[which(cdata$N==365 & cdata$Prop.miss<4),]


##Pull those sites and plot
low.miss.1<-subset(sp, site_name==low.miss$site_name[1] & year==2016)
low.miss.2<-subset(sp, site_name==low.miss$site_name[2] & year==2016)

##Simulate light for each site
hms<-rep("12:00:00", times=length(low.miss.1$date))
low.miss.1$date<-as.POSIXct(paste(low.miss.1$date.f, hms),tryFormats = c("%Y-%m-%d %H:%M:%OS"))
low.miss.2$date<-as.POSIXct(paste(low.miss.2$date.f, hms),tryFormats = c("%Y-%m-%d %H:%M:%OS"))

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
low.miss.1$sim.light<-lightest(time=low.miss.1$date, lat=sd$lat[sd$site_name==low.miss$site_name[1]], longobs=sd$lon[sd$site_name==low.miss$site_name[1]],longstd=75)
low.miss.2$sim.light<-lightest(time=low.miss.2$date, lat=sd$lat[sd$site_name==low.miss$site_name[2]], longobs=sd$lon[sd$site_name==low.miss$site_name[2]],longstd=75) 

##Impute missing covariates for prediction
#save location of missing data

lowmiss1_missing_integers<-which(is.na(low.miss.1$GPP))
lowmiss2_missing_integers<-which(is.na(low.miss.2$GPP))

#Do multiple imputations and save shortwave light and discharge  
lowmiss1<-as.data.frame(cbind(low.miss.1$GPP,low.miss.1$temp.water,low.miss.1$discharge, low.miss.1$shortwave,low.miss.1$date.f))
colnames(lowmiss1)<-c("GPP", "temp","q","meas.light", "ts")
z2<-amelia( lowmiss1, m = 5, p2s=1, ts="ts", lags="GPP", bounds=rbind(c(1,0,Inf)))
i1<-z2$imputations$imp1
i2<-z2$imputations$imp2
i3<-z2$imputations$imp3
i4<-z2$imputations$imp4
i5<-z2$imputations$imp5
arr <- abind(i1,i2,i3,i4,i5, along = 3)
z3<-as.data.frame(rowMeans(arr, dims = 2))
low.miss.1$imp.shortwave<-z3$meas.light
low.miss.1$imp.discharge<-z3$q



lowmiss2<-as.data.frame(cbind(low.miss.2$GPP,low.miss.2$temp.water,low.miss.2$discharge, low.miss.2$shortwave,low.miss.2$date.f, low.miss.2$sim.light))
colnames(lowmiss2)<-c("GPP", "temp","q","meas.light", "ts","sim.light")
z2<-amelia( lowmiss2, m = 5, p2s=1, ts="ts", lags="GPP", bounds=rbind(c(1,0,Inf)))
i1<-z2$imputations$imp1
i2<-z2$imputations$imp2
i3<-z2$imputations$imp3
i4<-z2$imputations$imp4
i5<-z2$imputations$imp5
arr <- abind(i1,i2,i3,i4,i5, along = 3)
z3<-as.data.frame(rowMeans(arr, dims = 2))
low.miss.2$imp.shortwave<-z3$meas.light
low.miss.2$imp.discharge<-z3$q
low.miss.2$imp.temp<-z3$temp


###Plot site data and light

ggplot(data=low.miss.1, aes(x=date, y=GPP))+
  geom_point(color="black", size=3)+
  geom_errorbar(aes(x=date,ymin=GPP.lower, ymax=GPP.upper), width=0.2, size=0.5)+
  geom_line(aes(x=date, y=imp.shortwave/100), color="orange")+
  geom_line(aes(x=date, y=imp.discharge/10), color="blue")+
  theme(legend.position="top")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  theme(legend.position="top")+
  ylab("GPP")+
  xlab("Time (days)")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggplot(data=low.miss.2, aes(x=date, y=GPP))+
  geom_point(color="black", size=3)+
  geom_errorbar(aes(x=date,ymin=GPP.lower, ymax=GPP.upper), width=0.2, size=0.5)+
  #geom_line(aes(x=date, y=imp.shortwave/100), color="orange")+
  geom_line(aes(x=date, y=(imp.discharge/3)-2), color="blue")+
  theme(legend.position="top")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  theme(legend.position="top")+
  ylab("GPP")+
  xlab("Time (days)")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())




# Force missing data ------------------------------------------------------


##Force some missing data-BY WEEK##
#Note-the model will still work if no data is missing but the objects in the next block of code (i.e., "y_index_mis", etc.)
#still need to be created
set.seed(50)

x<-low.miss.2$GPP
N<-length(x)
y_miss<-list(x,x,x,x,x,x,x,x)
#vector of missing number amounts
missing_n<-c(0,round(N*0.15),round(N*0.3),round(N*0.4),round(N*0.5),round(N*0.6),round(N*0.7),round(N*0.8))
#Create a list to store missing data integers
y_missing_integers<-c()
#Create a vector of the length of data representing the indexes
index<-c()
index[[1]]<-1:N

#Create initial missing data integer set. Then build on with the for loop
y_missing_integers[[1]]<-which(is.na(y_miss[[1]]))
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
  z[1]<-x[1]
  y_miss[[i]]<-z
}
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
for(i in 1:length(missing_n)){
  r<-y_miss[[i]]
  r[y_missing_integers[[i]]]<--100 #arbitrary number
  r[1]<-x[1]
  y_miss[[i]]<-r
} 

summary(y_miss)

#saveRDS(low.miss.1, file = "low.miss.1.missing.RDS") 
#saveRDS(low.miss.2, file = "low.miss.2.missing.RDS") 


# Stan model --------------------------------------------------------------


sink("toy2p.stan")

cat("
 
/*----------------------- Data --------------------------*/
  /* Data block: defines the objects that will be inputted as data */
  data {
    int N; // Length of state and observation time series
    vector[N] y_miss; // Observations
    real z0; // Initial state value
    vector[N] light;
    vector[N] dis;
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
    real b2;
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
    b2 ~ normal(0,5);
    // Distribution for the first state
    y[1] ~ normal(z0, sdp);
    // Distributions for all other states
    for(t in 2:N){
      y[t] ~ normal(b0+y[t-1]*phi+light[t]*b1+b2*dis[t], sdp);
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
 y_rep[t]=normal_rng(b0+y[t-1]*phi+light[t]*b1+b2*dis[t], sdp);
 }
 
  }
  
    "
    ,fill=TRUE)
sink()
closeAllConnections()


# Bayesian Parameter estimation -------------------------------------------



fit.miss1 <- vector("list",length(missing_n))
model<-"toy2p.stan"
model<-stan_model(model)
for(i in 1:length(missing_n)){
  ##Load data
  data <- list(   N = length(y_miss[[i]]),
                  y_nMiss = unlist(y_nMiss[[i]]),
                  y_index_mis = unlist(y_index_mis[[i]]),
                  y_miss= unlist(y_miss[[1]]),
                  light=log(low.miss.2$imp.shortwave),
                  dis=log(low.miss.2$imp.discharge),
                  #temp=log(low.miss.2$imp.temp),
                  z0=y_miss[[i]][1])
  
  ##Run Stan
  fit.miss1[[i]]<- rstan::sampling(object=model, data = data,  iter = 4000, chains = 4, control=list(max_treedepth = 15))
}

#saveRDS(fit.miss1, file = "full_bayes_lowmiss2_withqandtemp.RDS") 


##HMC diagnostics
rstan::check_hmc_diagnostics(fit.miss1[[8]])

##Traceplots
traceplot(fit.miss1[[6]], pars=c("phi", "b0","sdp", "b1"))

##PPC
fit_extract<-rstan::extract(fit.miss1[[6]])
yrep1 <- fit_extract$y_rep
samp100 <- sample(nrow(yrep1), 100)
ppc_dens_overlay(y_miss[[6]], yrep1[samp100, ]) + xlim(-2, 15)
ppc_stat(y, yrep1[samp100, ])


##Exctract parameters
fit_summary_pars_bayes <- vector("list",length(missing_n))

for (i in 1:length(missing_n)){
  fit_summary_pars_bayes[[i]]<-(summary(fit.miss1[[i]], pars=c("sdp","phi", "b1", "b0", "b2"), probs=c(0.025,.5,.975))$summary)
}
#check
fit_summary_pars_bayes[[1]]

##Unlist,cleanup, and add factors
fit_summary_bayes<-as.data.frame(do.call(rbind.data.frame, fit_summary_pars_bayes)) #Unlist
fit_summary_bayes$param<-rep(c("sdp","phi", "b1","b0", "b2"), times=length(missing_n) )#add parameter factor
fit_summary_bayes$prop.missing<-rep(prop.miss, times=rep(5,times=length(missing_n))) #add prop of missing data
row.names(fit_summary_bayes)<-1:(length(missing_n)*5) #remove row names
fit_summary_bayes$param<-as.factor(fit_summary_bayes$param) #make factor
fit_summary_bayes$prop.missing<-as.factor(fit_summary_bayes$prop.missing) #make factor
colnames(fit_summary_bayes)<-c("mean", "se-mean", "sd", "min", "med", "high","n_eff", "rhat", "param","prop.missing")
summary(fit_summary_bayes) #check it looks good
head(fit_summary_bayes)

##Save lowest missing data as "known" for comparison with more missing data
known.data<-fit_summary_bayes$mean[fit_summary_bayes$prop.missing==3]
known<-as.data.frame(known.data)

known.param<-c("sdp","phi", "b1","b0", "b2")
#known.missing<-c(5, 5, 5,5)
known<-as.data.frame(cbind(known.data, known.param))
known$known.data<-as.numeric(as.character(known$known.data))


saveRDS(fit.miss1, file = "full_lowmiss2_day_bayes_light_Q.RDS")
saveRDS(fit_summary_bayes, file = "summary_lowmiss2_day_bayes_light_Q.RDS")
saveRDS(known, file = "known_lowmiss2_day_bayes_light_Q.RDS")


##Plot parameters
ggplot(data=fit_summary_bayes, aes(x=mean, y=prop.missing ))+
  geom_point(aes( color=param, group=param),size=3,position=position_dodge(0.5))+
  theme_classic()+
  geom_errorbar(aes(xmin=min, xmax=high,color=param, group=param), width=0.2, size=0.5,position=position_dodge(0.5))+
  theme(legend.position="top")+
  ylab("Percent of Missing Data")+
  xlab("Parameter Estimate")+
  scale_color_manual(values=c("blue", "green", "darkgray", "black", "red"))+
  scale_x_continuous(limits=c(-2,1.1))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  geom_vline(xintercept = known.data,color=c("red", "black", "green", "blue", "darkgray"))


##Create object with estimated data
fit_summary_bayes_data<-(summary(fit.miss1[[6]], probs=c(0.025,.5,.975))$summary)
my_data <- as_tibble(fit_summary_bayes_data)
y_est_data<-my_data %>% slice(1:length(y_index_mis[[6]]))

##Create object with observed data
date<-1:length(y_miss[[6]])
y_obs_data<-as.data.frame(cbind(date[y_index_mis[[6]]],low.miss.2$GPP[y_index_mis[[6]]]))
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


raw_data<-as.data.frame(cbind(low.miss.1$GPP,date))

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



# Multiple Imputations ----------------------------------------------------



data.amelia<- vector("list",length(missing_n))
for(i in 1:length(missing_n)){
  y_miss1<-na_if(y_miss[[i]], -100)
  y_miss1[i][[1]]<-y_miss[[i]][1]
  w<-as.data.frame(cbind(y_miss1,low.miss.2$date.f, low.miss.2$discharge, 
                           low.miss.2$imp.shortwave ))
  colnames(w)<-c("y_miss1", "ts", "discharge", "imp.light")
  z2<-amelia(w, m = 5, p2s=1, ts="ts", lags="y_miss1", bounds=rbind(c(1,0,Inf)))
  summary(z2)
  data.amelia[[i]]<-z2
}


#Assign NAs to the missing data in each list
for(i in 8){#:length(missing_n)){
  z1<-data.amelia[[i]]$imputations$imp1$y_miss1
  z2<-data.amelia[[i]]$imputations$imp2$y_miss1
  z3<-data.amelia[[i]]$imputations$imp3$y_miss1
  z4<-data.amelia[[i]]$imputations$imp4$y_miss1
  z5<-data.amelia[[i]]$imputations$imp5$y_miss1
  q1<-unique(z1[y_missing_integers[[i]]])
  q2<-unique(z2[y_missing_integers[[i]]])
  q3<-unique(z3[y_missing_integers[[i]]])
  q4<-unique(z4[y_missing_integers[[i]]])
  q5<-unique(z5[y_missing_integers[[i]]])
  q.ts<-unique(low.miss.2$date.f[y_missing_integers[[i]]])
  q<-cbind(q1,q2,q3,q4,q5)
  q.mean<- apply(q[,-1], 1, mean)
  q.min<-apply(q[,-1], 1, min)
  q.max<-apply(q[,-1], 1, max)
  q.sd<- apply(q[,-1], 1, sd)
  x_obs_data<-as.data.frame(cbind(x,low.miss.2$date.f))
  x_imp_data<-as_tibble(cbind(q.mean,q.sd,q.ts, q.min,q.max))
  x_imp_data<-as.data.frame(x_imp_data %>% slice(1:length(unique(y_missing_integers[[i]]))))
  x_obs<-as.data.frame(low.miss.2$GPP[unique(y_missing_integers[[i]])])
  x_combined<-as.data.frame(cbind(x_imp_data$q.mean,x_imp_data$q.sd,x_obs,x_imp_data$q.min,x_imp_data$q.max))
  colnames(x_combined)<-c("imp.mean", "imp.sd", "obs", "imp.min","imp.max")
  g<-ggplot(data=x_imp_data, aes(x=q.ts, y=q.mean))+
    geom_errorbar(data=x_imp_data,aes(x=q.ts,ymin=q.min, ymax=q.max))+
    geom_point(data=x_obs_data, aes(x=V2, y=x, color="gray"), size=3)+
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
  
  s<-ggplot(data=x_combined, aes(x=obs, y=imp.mean))+
    geom_abline(intercept=0, slope=1)+
    theme_classic()+
    geom_errorbar(aes(ymin=imp.min, ymax=imp.max))+
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





fit.amelia <- vector("list",length(missing_n))
data.stan.amelia <- vector("list",(length(missing_n)))
list.1<- vector("list",(5))
for(i in 1:length(missing_n)){##two because the list holding the lists starts with a list of misc data that arent the estimates
  for(g in 1:5){
    list.1[[g]] <- data.amelia[[i]]$imputations[[g]]$y_miss1
    
  }
  data.stan.amelia[[i]]<-list.1
}

plot(data.stan.amelia[[5]][[2]],data.stan.amelia[[5]][[3]])
##Run Stan


fit.stan.miss.amelia2 <- vector("list",length(missing_n))
fit.stan.miss.imp<- vector("list",5)
model<-"toy2p.stan"
model<-stan_model(model)
for(i in 1:3){#:length(missing_n)){
  for(g in 1:5){
    ##Load data
    data <- list(   N = length(y_miss[[i]]),
                    y_nMiss = unlist(y_nMiss[[1]]),
                    y_index_mis = unlist(y_index_mis[[1]]),
                    y_miss= unlist(data.stan.amelia[[i]][[g]]),
                    light=log(low.miss.2$imp.shortwave),
                    dis=log(low.miss.2$imp.discharge),
                    z0=y_miss[[i]][1])
    
    ##Run Stan
    fit.stan.miss.imp[[g]]<- rstan::sampling(object=model, data = data,  iter = 4000, chains = 4, control=list(max_treedepth = 15))
  }
  fit.stan.miss.amelia2[[i]]<-fit.stan.miss.imp
}

##Pull param estimates into list
fit_summary_pars_amelia2 <- vector("list",length(missing_n))
list.2<- vector("list",5)
for (i in 1:length(missing_n)){
  for (g in 1:5){
    list.2[[g]] <-summary(fit.stan.miss.amelia2[[i]][[g]], pars=c("sdp","phi", "b1", "b0", "b2"), probs=c(0.025,.5,.975))$summary
  }
  fit_summary_pars_amelia2[[i]]<-list.2
}
fit_summary_pars_amelia2[[2]]#<-fit_summary_pars_amelia[-1]

##Unlist,cleanup, and add factors
fit_summary<- vector("list",length(missing_n))

for(i in 1:length(missing_n)){
  fit_summary[[i]]<-as.data.frame(do.call(rbind,fit_summary_pars_amelia2[[i]]))
  fit_summary[[i]]$param<-rep(c("sdp","phi", "b1","b0","b2"), times=5 )#add parameter factor
  fit_summary[[i]]$prop.missing<-rep(prop.miss[i], each=5*5) #add prop of missing data
  fit_summary[[i]]$imp.set<-rep(1:5, each=5) #add prop of missing data
  row.names(fit_summary[[i]])<-1:25 #remove row names
  fit_summary[[i]]$param<-as.factor(fit_summary[[i]]$param) #make factor
  fit_summary[[i]]$prop.missing<-as.factor(fit_summary[[i]]$prop.missing) #make factor
  colnames(fit_summary[[i]])<-c("mean", "se-mean", "sd", "min", "med", "high","n_eff", "rhat", "param","prop.missing", "imp.set")
}
fit_summary<-as.data.frame(do.call(rbind,fit_summary))


summary(fit_summary) #check it looks good
head(fit_summary)

summ<-ddply(fit_summary, c("param", "prop.missing"), summarize,
            avg    = mean(mean),
            ci.min = mean(min),
            ci.max = mean(high))

#,vw= sum(se.mean.sq)/5)

fit_summary<-left_join(fit_summary, summ, by=c("param", "prop.missing"))


##Save lowest missing data as "known" for comparison with more missing data
known.data<-subset(fit_summary,prop.missing==3 & imp.set==1)

known.param<-c("sdp","phi", "b1","b0", "b2")
#known.missing<-c(5, 5, 5,5)
known<-as.data.frame(cbind(known.data$mean, known.param))
colnames(known)<-c("mean", "param")
known$mean<-as.numeric(as.character(known$mean))

ggplot(data=fit_summary, aes(x=avg, y=prop.missing ))+
  geom_point(aes( color=param, group=param),size=3,position=position_dodge(0.5))+
  theme_classic()+
  geom_errorbar(aes(xmin=ci.min, xmax=ci.max,color=param, group=param), width=0.2,position=position_dodge(0.5))+
  theme(legend.position="right")+
  ylab("Percent of Missing Data")+
  xlab("Parameter Estimate")+
  scale_color_manual(values=c("blue", "green",  "red","darkgray", "black"))+
  scale_x_continuous(limits=c(-5,1.5), breaks= c(-5,-4,-3,-2,-1,0,1,1))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  geom_vline(xintercept = known$mean,color=c("darkgray", "black", "green", "blue", "red"))

pairs(fit.stan.miss.amelia2[[8]][[5]], pars=c("phi", "b0","sdp", "b1"))


saveRDS(fit.stan.miss.amelia2, file = "full_amelia_lowmiss2_model_var.RDS")
saveRDS(known, file = "known_data_amelia_lowmiss2_model_var.RDS") 
saveRDS(fit_summary, file = "summary_amelia_lowmiss2_model_var.RDS") 


filename <- file.choose()
Canteen_clean <- readRDS(filename)



# GAM ---------------------------------------------------------------------


library(mgcv)
library(MASS)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(feather)

x <- seq(0, pi * 2, 0.1)
sin_x <- sin(x)
y <- sin_x + rnorm(n = length(x), mean = 0, sd = sd(sin_x / 2))
Sample_data <- data.frame(y,x)

ggplot(Sample_data, aes(x, y)) + geom_point()

gam_y <- gam(y ~ s(x), method = "REML")
x_new <- seq(0, max(x), length.out = 100)
y_pred <- predict(gam_y, data.frame(x = x_new))
ggplot(Sample_data, aes(x, y)) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x))

par(mfrow = c(2,2))
gam.check(gam_y)

low.miss.1$doy<-seq(1:length(low.miss.1$GPP))

gam_1 <- gam(GPP ~ s(doy, bs = "cc", k = 5)+
               s(sim.light,bs = "re", k = 5),
             data = low.miss.1,
             family = gaussian,
             method = "REML")

plot(gam_1, shade = TRUE)


summary(gam_1)
