# Load packages
library(here)
library(tidyverse)
library(rstan)
library(Amelia)
library(forecast)
library(xts)
library(nlme)
library(tidyverse)
library(lubridate)


# Reference A. Stears' code with helpful function for removing data
source("functions/missing_data_functions.R")


########## Run code for ARIMA methods over GPP 1 simulated dataset w/ increasing missingness
## 3 diffnet types (MAR low auto, MAR high auto, and MNAR) #####

## ARIMA method 1: Drop NAs ###

source("Functions/Arima_drop_function.R")

source("Model Runs/Data_Deletion_ARIMA.R")

head(paramdropall)

## ARIMA method 2: Kalman filter (preserves the NAS) ###

source("Functions/Arima_Kalman_function.R")

source("Model Runs/Kalman_ARIMA.R")

head(paramKalmanall)

## ARIMA method 3: Multiple imputations with AMELIA ###


source("Functions/Arima_MI_function.R")

source("Model Runs/MI_ARIMA.R")

head(paramMIall)

########## Run code for ARIMA methods over GPP 1 simulated dataset w/ increasing missingness #####

source("Functions/model_fitting.R")

source("Model Runs/DA_BRMS_GPP.R") ## COMPUTER CRASHES WHEN i RUN THIS ALL AT ONCE ##

##ran them seperately and saved .csv. combining them all together back in here #

BRMSall<-rbind(BRMSMARlowauto, BRMSMARhighauto, BRMSMNAR)

head(BRMSall)


##### FORMATTING ALL ###

### Estimates without missing data ####

gauss_sim_MAR_datasets <- readRDS("data/Missingdatasets/gauss_sim_randMiss_autoCorr_10.rds")

## no missing###

sim1<-gauss_sim_MAR_datasets [[1]][["y"]][["y_noMiss"]]

realphi<-gauss_sim_MAR_datasets[[1]][["sim_params"]][["phi"]]
realbetas<-gauss_sim_MAR_datasets[[1]][["sim_params"]][["beta"]]

trueestdf <- data.frame (param = c("ar1", "intercept", "light", "discharge"), value = c(realphi, realbetas))

covariates<-gauss_sim_MAR_datasets[[1]][["sim_params"]][["X"]]

covariatesX<-as.matrix(covariates[,2:3])

days<-seq(1, 365)

sim1df<-as.data.frame(cbind(days=days, GPP=sim1, light=covariates[,2], discharge=covariates[,3]))


##quick plot simulated dataset##
simGPP <- ggplot(sim1df, aes(x=days, y=GPP))+
  geom_point(size = 2, color="chartreuse4") + 
  geom_line(size = 1, color="chartreuse4")+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title="Simulated GPP")


#########################################################################################################################
## ARIMA estimate with no missing data ##
###########################################################################################################################

##matrix of covariates##
X = matrix(c(sim1df$light,  sim1df$discharge), ncol = 2)

Arimanomissing=arima(sim1df$GPP, order = c(1,0,0), xreg = X)

Arimanomissing

arimaestdf <- data.frame (param = c("ar1", "intercept", "light", "discharge"), value =Arimanomissing$coef)



###################  PLOTTING #######################

#### plotting ARIMA AND DATA AUGMENTATION (STAN) TOGETHER ###


##########arima output#############

#combine all the arima outputs ##

paramallARIMA<-rbind(paramdropall,paramKalmanall,paramMIall)

paramallARIMA<-paramallARIMA %>% select(type,missingness, param, value, SE, missingprop,autocorr)

## change AR 1 to be phi ##

paramallARIMA2<-paramallARIMA

paramallARIMA2$param <- str_replace(paramallARIMA2$param, "ar1", "phi")

## distinguish between arima and stan ###

paramallARIMA2$type <- str_replace(paramallARIMA2$type, "Data Deletion", "Data Deletion: Arima")
paramallARIMA2$type <- str_replace(paramallARIMA2$type, "Multiple imputations", "Multiple imputations: Arima")
paramallARIMA2$type <- str_replace(paramallARIMA2$type, "Kalman filter", "Kalman filter: Arima (default)")

#### add ARIMA and STAN dataframes together ###

paramallARIMASTAN<-rbind(paramallARIMA2, BRMSall[,2:8])

## drop sdp for comparison with arima ###
paramallARIMASTAN2<- paramallARIMASTAN %>% filter(param != "sigma" ) 



#paramallARIMASTAN2$type <- factor(paramallARIMASTAN2$type, levels=c("Data Deletion: Arima", "Kalman filter: Arima (default)", "Multiple imputations: Arima", "Data augmentation: STAN"), ordered=TRUE)

#arimastanMAR<-ggplot(data=paramallARIMASTAN2, aes(x=as.numeric(missingprop), y=value))+
  #facet_grid(~factor(param, levels=c("intercept", "phi", "light", "discharge"),exclude = NA)~ type, scales="free_y")+
  #geom_hline(data=trueestdf2, aes(yintercept=value), colour="gray")+
  #geom_hline(data=arimaestdf, aes(yintercept=value), colour="light blue")+
  #geom_point(size=1)+
  #geom_errorbar(aes(ymin=value-SE, ymax=value+SE), size=0.3)+
  #theme_bw()+
  #xlab("Percent of Missing Data (Missing at Random)")+
  #ylab("Parameter estimate")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 6))


#pdf(file = "Figures/MARcomparison.pdf",width = 8, height = 5)
#arimastanMAR
#dev.off()

trueestdf2 <- data.frame (param = c("phi", "intercept", "light", "discharge"), value = c(realphi, realbetas))

paramallARIMASTAN2$missingness <- factor(paramallARIMASTAN2$missingness, levels = c("MAR low auto", "MAR high auto", "MNAR"))

# filter missingness > 0.77 ##

paramallARIMASTAN3<-paramallARIMASTAN2 %>% filter(missingprop < 0.77)%>% filter(missingprop!= 0.22)%>% filter(missingprop!= 0.24)



### different version of the same figure ##

allmethodsGPP<-ggplot(data=paramallARIMASTAN3, aes(x=as.numeric(missingprop), y=value, color=type))+
  facet_grid(~factor(param, levels=c("intercept", "phi", "light", "discharge"),exclude = NA)~ missingness, scales="free_y")+
  geom_hline(data=trueestdf2, aes(yintercept=value), colour="gray")+
  #geom_hline(data=arimaestdf, aes(yintercept=value), colour="light blue")+
  geom_point(size=0.75, position = position_dodge(width=0.03))+
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), size=0.3, width=0, position = position_dodge(width=0.03))+
  theme_bw()+
  xlab("Proportion of missing data")+ theme(legend.position="top")+theme(legend.title=element_blank())+
  ylab("Parameter estimate")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))


pdf(file = "Figures/MARautoMNAR.pdf",width = 8, height = 5)
allmethodsGPP
dev.off()




####

