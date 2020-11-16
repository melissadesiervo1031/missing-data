##Load library
library(Amelia)
library(readr)
library(tidyr)
library(ggplot2)
library(chron)
library(scales)
library(dplyr)
library(data.table)
library(lubridate)
##Load tutorial data
data(freetrade)

##Summarize tutorial data. Note variables with missing data
summary(freetrade)

###Impute missing data. 
##STANDARD
#'m'= number of imputed datasets
#'#'p2s'=0 prints nothing, 1 prints some output, 2 prints a lot of stuff including chain status
#'idvars'=VARIABLE(s). Variables to be retained but are not included in the imputation

##DATA DESCRIPTION
#'ords'=VARIABLE will turn all variables into oridinal
#'noms'=VARIABLE. Makes sure that categoircal dummy variables are imputed only as categories that exist in the data
#'logs'=VARIABLE. Natural log transformation
#'sqrts'=VAIABLE. Square root transformation
#'lgstc'=VARIABLE. Bounds data for proportions between 0 and 1

##TIME OR SPATIAL CONSIDERATIONS
#'ts'= time series variable. Useless without command also for polytime and intercs
#'polytime'= dictates the polynomial order (1-3)
#'lags'=VARIABLE. includes a lag of the variable in the imputation
#'leads'=VARIABLE. includes a lead of the variable in the imputation
#'cs'= spatial data variable. Useless without command also for polytime and intercs
#'intercs' = TRUE. Allows polynomial fitting to happen for each cross section factor
#'tscsPlot'. FUNCTION to plot time series with observed and imputed data

###PRIORS
##RIDGE PRIORS
#When there is a high degree of missingmess, very strong correlation among variables, or
#when the observations is slightly greater than the number of parameters (p(p+3)/2; p= variables)
#the analysis model will depend on the choice of imputation model
#Adds stability for the sake of some accuracy.
#'emperi'=VALUE OF OBSERVATIONS. 0.5-1% of observations is normal, 5% is moderate, 10% upper limit.(e.g., .01*nrow(VARIABLE); for 1%)

##OBSERVATION PRIORS
#Priors on individual data points based on previous knowedge or personal experience.
#On individual obervation data points NOT the model parameters

##Shatto Data
setwd("~/Data/SHAT/metabolism modeling/R_code_for modeling")

##Load data 
#metabolism estimates
alldata<-read.csv("All Metabolism_Metabolizer_2008-2017_61919_CTL.csv")
#format date
alldata$Date<-as.Date(alldata$Date, format="%m/%d/%Y")
#GPP greater than 0
alldata<-subset(alldata,GPP >= 0)

#fill in missing dates
alldata_complete<-complete(Date = seq.Date(min(alldata$Date), max(alldata$Date), by="day"), data=alldata)

#Discharge
mydataQ<-read.csv("ShattoQ.csv")
mydataQ$Date<-as.Date(mydataQ$Date,format="%m/%d/%Y")

#Light
mydata.light<-read.csv("SDW_estimated_light.csv")
mydata.light$date<-as.Date(mydata.light$date)
colnames(mydata.light)<-c("Date", "light")

#Temp
mydata.temp<-read.csv("SDW.temp.csv")
mydata.temp$Date<-as.Date(mydata.temp$Date,format="%m/%d/%Y")

#Merge Discharge with all data
setDT(mydataQ)
alldata_complete_Q <- left_join(alldata_complete, mydataQ, by="Date")

#Merge light with all data
setDT(mydata.light)
alldata_complete_Q_L <- left_join(alldata_complete_Q, mydata.light, by="Date")

#Merge temp with all data
setDT(mydata.temp)
alldata_complete_Q_L_T <- left_join(alldata_complete_Q_L, mydata.temp, by="Date")

#Retain important columns
alldata_complete_Q_L_T<-alldata_complete_Q_L_T[,c(1,6,8,19,20,21)]



##Add seasons based on month
seasons = function(x){
  if(x %in% 4:6) return("Spring")
  if(x %in% 7:9) return("Summer")
  if(x %in% 10:12) return("Fall")
  if(x %in% c(1:3)) return("Winter")
  }

alldata_complete_Q_L_T$Season=sapply(month(alldata_complete_Q_L_T$Date), seasons)
#alldata_complete_Q_L$Season=as.factor(alldata_complete_Q_L$Season)

##Add Year
alldata_complete_Q_L_T$Year=year(alldata_complete_Q_L_T$Date)

##Make a dataframe
alldata_complete_Q_L_T<-as.data.frame(alldata_complete_Q_L_T)

##Summarize with the number of NA's
summary(alldata_complete_Q_L_T)

##Basic Impute
a.out <- amelia(alldata_complete_Q_L_T, m = 5, p2s=1,idvars=c("Date","Season"))
a.out

##Plot imputed vs observed
plot(a.out,which.vars=3)


##TS Impute
a.out.ts <- amelia(alldata_complete_Q_L_T, m = 5, p2s=1,cs="Year", ts="Date", polytime = 1, intercs = FALSE,idvars=c("Season"),  emperi=.1*nrow(Date))#
a.out.ts

tscsPlot(a.out.ts, var=3, cs="Season", ts="Date",plotall=TRUE)

##Plot imputed vs observed
plot(a.out.ts,which.vars=3)

write.amelia(obj=a.out.ts, file.stem = "outdata.ts", separate=FALSE, orig.data=TRUE)
sapply(a.out.ts$imputations[[3]],  function(y) sum(is.na(y)))








#Stream metabolism Example
#From Appling et al. 2018 dataset
setwd("C:/Users/mtrentman/IDrive-Sync/Postdoc/Estimating missing data/daily_predictions")

#Site meta data
site<-read_tsv("site_data.tsv")

#metabolism data
daily_prediction<-read_tsv("daily_predictions.tsv")
daily_prediction$site_name<-as.factor(daily_prediction$site_name)

summary(daily_prediction)

#filter one site--Devils Tower
daily_devil<-subset(daily_prediction, site_name=="nwis_02266300")

#Retain only important variables for Amelia (more data will slow it way down)
sub_daily<-daily_devil[,c(1,3,4,9,14,24,26,27)]

##Fill in dates for missing data
##Fill in missing dates
sub_daily_complete<-complete(date = seq.Date(min(sub_daily$date), max(sub_daily$date), by="day"), data=sub_daily)

##Summarize with the number of NA's
summary(sub_daily_complete)

a.out <- amelia(sub_daily_complete, m = 5, p2s=1,idvars=c("site_name", "date"))
a.out

##Histogram of the i imputed data (3 here)
hist(a.out$imputations[[3]]$tariff, col="grey", border="white")

##Save the imputations to an R file
save(a.out, file = "imputations.RData")

##Write a CSV for each imputed dataset
write.amelia(obj=a.out, file.stem = "outdata")

##Write a stacked CSV of imputed datasets
#Creates impvar column to indicate the dataset # of imputations (i.e., for m=5, values will be 1-5)
write.amelia(obj=a.out, file.stem = "outdata", separate=FALSE, orig.data=TRUE)
