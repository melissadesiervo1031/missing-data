##Load library
library(Amelia)
library(readr)
library(tidyr)
library(ggplot2)
library(chron)
library(scales)
library(dplyr)
library(data.table)

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
