library(tidyr)
library(MASS)
#Time series
setwd("G:/My Drive/Postdoc/Metabolism time series")

#Load Data
ts.SDW<-read.csv("SDW.ts_with light.csv")
ts.SDW$Year<-as.factor(ts.SDW$year)

##Format date for filling in dates with missing data
ts.SDW$date<-as.Date(ts.SDW$date,format="%m/%d/%Y")

##Remove negative GPP values
ts.SDW<- subset(ts.SDW,GPP >= 0)

##Fill in missing dates
df<-complete(date = seq.Date(min(ts.SDW$date), max(ts.SDW$date), by="day"), data=ts.SDW)

##Code data to be missing=1 or observed=0
N<-ifelse(is.na(df$GPP),1,0)

## Calculate length of data gaps
#prepare vectors for loop
x<-NA 
x[1]<-N[1]
length_miss[1]<-1

#loop over data and return the length of gaps
for (i in 2:length(df$GPP)){
x[i]<-ifelse(N[i]==1,1+x[i-1], 0) #If data is missing add to yesterdays value by 1. If data observed then 0.
length_miss[i]<-ifelse(x[i]==0,x[i-1], NA) #Save length of gap on the day data is observed (0) after a gap
length_miss[i]<-ifelse(length_miss[i]==0,NA,length_miss[i]) #Replace days with observed data with NA
}

## Calculate length of observed data before gap. Same as above but for lengths of observed data
z<-NA
z[1]<-N[1]
length_obs<-NA
for (i in 2:length(df$GPP)){
  z[i]<-ifelse(N[i]==0,1+z[i-1], 0)
 length_obs[i]<-ifelse(z[i]==0,z[i-1], NA)
  length_obs[i]<-ifelse(length_obs[i]==0,NA,length_obs[i] )
}


## save to full dataset and check that everything is correct
df$length_obs<-length_obs
df$length_miss<-length_miss

## estimate negative binomial distribution parameters of missing data length integers
length_miss_parm<-length_miss[complete.cases(length_miss)]
ff <- fitdistr(length_miss_parm, "Negative Binomial")
ff

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


