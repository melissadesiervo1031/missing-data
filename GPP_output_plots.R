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
# makeMissing()
source("missing_data_functions.R")




############# MISSING COMPLETELY AT RANDOM (MCAR)#########################################


#### Read in the missing dataframes that Alice S. made #####
###### Pull out the first in 1000 nested lists for this code ### (Eventually loop over all the lists)

gauss_sim_MAR_datasets <- readRDS(here("data/Missingdatasets/gauss_sim_randMiss.rds"))


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



##For nested list of GPP datasets with increasing MAR data add back in the date column and the covariates## 

GPP_sim_MAR<- gauss_sim_MAR_datasets [[1]][["y"]]

GPP_sim_MAR_2 <-lapply(X = GPP_sim_MAR, FUN = function(X)   cbind.data.frame(GPP=X, days=sim1df$days, light = sim1df$light, discharge = sim1df$discharge))



#########################################################################################################################
## ARIMA estimate with no missing data ##
###########################################################################################################################

##matrix of covariates##
X = matrix(c(sim1df$light,  sim1df$discharge), ncol = 2)

Arimanomissing=arima(sim1df$GPP, order = c(1,0,0), xreg = X)

Arimanomissing

arimaestdf <- data.frame (param = c("ar1", "intercept", "light", "discharge"), value =Arimanomissing$coef)

#### PLOT ALL THE ARIMA METHODS TOGETHER ###############


paramdroplong2 
paramMIlong2
paramNAlong2 

paramallARIMA<-rbind(paramdroplong2, paramMIlong2, paramNAlong2)


arimallmethods<-ggplot(data=paramallARIMA, aes(x=as.numeric(missingprop), y=value))+
  facet_grid(~factor(param, levels=c("ar1", "intercept", "light", "discharge"))~ type, scales="free")+
  geom_hline(data=trueestdf, aes(yintercept=value), colour="salmon")+
  geom_hline(data=arimaestdf, aes(yintercept=value), colour="light blue")+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE))+
  theme_bw()+
  ggtitle("Parameter estimates generated from Arima models")+
  xlab("Percent of Missing Data (Missing at Random)")+
  ylab("Parameter estimate")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#### plotting ARIMA AND DATA AUGMENTATION (STAN) TOGETHER ###

## change AR 1 to be phi ##

paramallARIMA2<-paramallARIMA

paramallARIMA2$param <- str_replace(paramallARIMA2$param, "ar1", "phi")

## distinguish between arima and stan ###

paramallARIMA2$type <- str_replace(paramallARIMA2$type, "Data Deletion", "Data Deletion: Arima")
paramallARIMA2$type <- str_replace(paramallARIMA2$type, "Multiple imputations", "Multiple imputations: Arima")
paramallARIMA2$type <- str_replace(paramallARIMA2$type, "Kalman filter", "Kalman filter: Arima (default)")

#### add ARIMA and STAN dataframes together ###

## ERROR (SE for ARIMA, SD for STAN) ###


paramallARIMA3 <- paramallARIMA2 %>% rename("error" = "SE") %>% mutate(errortype="SE")
DAparamdf4 <- DAparamdf3 %>% rename("error" = "SD") %>% mutate(errortype="SD")


paramallARIMASTAN<-rbind(paramallARIMA3, DAparamdf4)

## drop sdp for comparison with arima ###
paramallARIMASTAN2<- paramallARIMASTAN %>% filter(param != "sdp" ) 

trueestdf2 <- data.frame (param = c("phi", "intercept", "light", "discharge"), value = c(realphi, realbetas))


#paramallARIMASTAN2$type <- factor(paramallARIMASTAN2$type, levels=c("Data Deletion: Arima", "Kalman filter: Arima (default)", "Multiple imputations: Arima", "Data augmentation: STAN"), ordered=TRUE)

arimastanMAR<-ggplot(data=paramallARIMASTAN2, aes(x=as.numeric(missingprop), y=value, color=errortype))+
  facet_grid(~factor(param, levels=c("intercept", "phi", "light", "discharge"),exclude = NA)~ type, scales="free_y")+
  geom_hline(data=trueestdf2, aes(yintercept=value), colour="gray")+
  #geom_hline(data=arimaestdf, aes(yintercept=value), colour="light blue")+
  geom_point(size=1)+
  geom_errorbar(aes(ymin=value-error, ymax=value+error), size=0.3)+
  theme_bw()+
  xlab("Percent of Missing Data (Missing at Random)")+
  ylab("Parameter estimate")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 6))


pdf(file = "Figures/MARcomparison.pdf",width = 8, height = 5)
arimastanMAR
dev.off()


################################ MISSING NOT AT RANDOM (MNAR) ##########################

## no missing###

sim1MNAR<-gauss_sim_MNAR_datasets [[1]][["y"]][["y_noMiss"]]

realphiMNAR<-gauss_sim_MNAR_datasets[[1]][["sim_params"]][["phi"]]
realbetasMNAR<-gauss_sim_MNAR_datasets[[1]][["sim_params"]][["beta"]]

trueestdfMNAR <- data.frame (param = c("ar1", "intercept", "light", "discharge"), value = c(realphiMNAR, realbetasMNAR))

covariates<-gauss_sim_MNAR_datasets[[1]][["sim_params"]][["X"]]

covariatesX<-as.matrix(covariates[,2:3])

days<-seq(1, 365)

sim1MNARdf<-as.data.frame(cbind(days=days, GPP=sim1MNAR, light=covariates[,2], discharge=covariates[,3]))

#########################################################################################################################
## ARIMA estimate with no missing data ##
###########################################################################################################################

##matrix of covariates##
X = matrix(c(sim1df$light,  sim1df$discharge), ncol = 2)

Arimanomissing=arima(sim1df$GPP, order = c(1,0,0), xreg = X)

Arimanomissing

arimaestdf <- data.frame (param = c("ar1", "intercept", "light", "discharge"), value =Arimanomissing$coef)


#### PLOT ALL THE ARIMA METHODS TOGETHER ###############


paramdroplongMNAR2 
paramMIlongMNAR2
paramNAlongMNAR2 

paramallARIMAMNAR<-rbind(paramdroplongMNAR2, paramMIlongMNAR2, paramNAlongMNAR2)


arimallmethodsMNAR<-ggplot(data=paramallARIMAMNAR, aes(x=as.numeric(missingprop), y=value))+
  facet_grid(~factor(param, levels=c("ar1", "intercept", "light", "discharge"))~ type, scales="free")+
  geom_hline(data=trueestdf, aes(yintercept=value), colour="salmon")+
  geom_hline(data=arimaestdf, aes(yintercept=value), colour="light blue")+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE))+
  theme_bw()+
  ggtitle("Parameter estimates generated from Arima models")+
  xlab("Percent of Missing Data (Missing NOT at Random)")+
  ylab("Parameter estimate")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#### plotting ARIMA AND DATA AUGMENTATION (STAN) TOGETHER ###

## change AR 1 to be phi ##

paramallARIMAMNAR2<-paramallARIMAMNAR

paramallARIMAMNAR2$param <- str_replace(paramallARIMAMNAR2$param, "ar1", "phi")

## distinguish between arima and stan ###

paramallARIMAMNAR2$type <- str_replace(paramallARIMAMNAR2$type, "Data Deletion", "Data Deletion: Arima")
paramallARIMAMNAR2$type <- str_replace(paramallARIMAMNAR2$type, "Multiple imputations", "Multiple imputations: Arima")
paramallARIMAMNAR2$type <- str_replace(paramallARIMAMNAR2$type, "Kalman filter", "Kalman filter: Arima (default)")

#### add ARIMA and STAN dataframes together ###

paramallARIMASTANMNAR<-rbind(paramallARIMAMNAR2, DAparamdfMNAR3)

## drop sdp for comparison with arima ###
paramallARIMASTANMNAR2<- paramallARIMASTANMNAR %>% filter(param != "sdp" ) 

trueestdf2 <- data.frame (param = c("phi", "intercept", "light", "discharge"), value = c(realphi, realbetas))

arimastanMAR

paramallARIMASTANMNAR2$type <- factor(paramallARIMASTANMNAR2$type, levels=c("Data Deletion: Arima", "Kalman filter: Arima (default)", "Multiple imputations: Arima", "Data augmentation: STAN"), ordered=TRUE)

arimastanMNAR<-ggplot(data=paramallARIMASTANMNAR2, aes(x=as.numeric(missingprop), y=value))+
  facet_grid(~factor(param, levels=c("intercept", "phi", "light", "discharge"),exclude = NA)~ type, scales="free_y")+
  geom_hline(data=trueestdf2, aes(yintercept=value), colour="salmon")+
  #geom_hline(data=arimaestdf, aes(yintercept=value), colour="light blue")+
  geom_point(size=1)+
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), size=0.3)+
  theme_bw()+
  xlab("Percent of Missing Data (Missing NOT at Random)")+
  ylab("Parameter estimate")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 6))


pdf(file = "Figures/MNARcomparison.pdf",width = 8, height = 5)
arimastanMNAR
dev.off()





