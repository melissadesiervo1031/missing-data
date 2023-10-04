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


# read in data resulting from beartooth
gauss_real_ModelResults <- readRDS("./data/model_results/gauss_real_ModelResults.rds")
gauss_sim_ModelResults <- readRDS("./data/model_results/gauss_sim_ModelResults.rds")

# make columns for "autocor" and "missingness"
gauss_real_ModelResults$autoCor <- gauss_real_ModelResults$missingprop_autocor %>% 
  str_extract(pattern = "0.[0-9]+$") %>% 
  as.numeric()
gauss_real_ModelResults[gauss_real_ModelResults$missingness=="MNAR", "autoCor"] <- NA
gauss_real_ModelResults$amtMiss <- gauss_real_ModelResults$missingprop_autocor %>% 
  str_extract(pattern = "0.[0-9]+") %>% 
  as.numeric

gauss_sim_ModelResults$autoCor <- gauss_sim_ModelResults$missingprop_autocor %>% 
  str_extract(pattern = "0.[0-9]+$") %>% 
  as.numeric()
gauss_sim_ModelResults[gauss_sim_ModelResults$missingness=="MNAR", "autoCor"] <- NA
gauss_sim_ModelResults$amtMiss <- gauss_sim_ModelResults$missingprop_autocor %>% 
  str_extract(pattern = "0.[0-9]+") %>% 
  as.numeric
gauss_sim_ModelResults <- gauss_sim_ModelResults %>% 
  mutate(value = as.numeric(value),
         SE = as.numeric(SE))

# Gaussian simulated data -------------------------------------------------
# get just "high autocorrelation MAR", "low autocorrelation MAR", and "MNAR" (for all simulations)
gauss_sim_figDat <- gauss_sim_ModelResults %>% 
  filter(missingness == "MNAR" # get data that is MNAR 
         | (missingness == "MAR" & autoCor <= 0.1) # get low autocor MAR data
           | (missingness == "MAR" & autoCor >= 0.9 )) 

gauss_sim_figDat[gauss_sim_figDat$autoCor <=0.1 & !is.na(gauss_sim_figDat$autoCor), "missingness"] <- "MAR_lowAutoCor"
gauss_sim_figDat[gauss_sim_figDat$autoCor  >= 0.9 & !is.na(gauss_sim_figDat$autoCor), "missingness"] <- "MAR_highAutoCor"

# calculate the standardized difference between parameter estimates and simulated values??
# for phi

gauss_sim_figDat <- gauss_sim_figDat %>% 
  dplyr::mutate(paramDiff_phi = (value - phi_sim)/phi_sim,
                paramDiff_intercept = (value - intercept_sim)/intercept_sim, 
                paramDiff_light = (value-light_sim)/light_sim,
                paramDiff_discharge = (value -discharge_sim)/discharge_sim)
  
# remove unnecessary values
gauss_sim_figDat[gauss_sim_figDat$param == "sigma", c("paramDiff_phi", "paramDiff_intercept", "paramDiff_light", "paramDiff_discharge")] <- NA
gauss_sim_figDat[gauss_sim_figDat$param == "phi", c("paramDiff_intercept", "paramDiff_light", "paramDiff_discharge")] <- NA
gauss_sim_figDat[gauss_sim_figDat$param == "light", c("paramDiff_phi", "paramDiff_intercept", "paramDiff_discharge")] <- NA
gauss_sim_figDat[gauss_sim_figDat$param == "intercept", c("paramDiff_phi", "paramDiff_light", "paramDiff_discharge")] <- NA
gauss_sim_figDat[gauss_sim_figDat$param == "discharge", c("paramDiff_phi", "paramDiff_intercept", "paramDiff_light")] <- NA

gauss_sim_figDat$paramDiff_stnd <- apply(X = gauss_sim_figDat[,c("paramDiff_phi", "paramDiff_intercept", "paramDiff_light", "paramDiff_discharge")],
                                        MARGIN = 1, 
                                        FUN = function(x) sum(x, na.rm = TRUE))
gauss_sim_figDat <- gauss_sim_figDat %>% 
  select(-paramDiff_phi,
         -paramDiff_intercept,
         -paramDiff_light,
         -paramDiff_discharge)
# make bins of proportion missingness
gauss_sim_figDat$amtMiss_bin <- round(gauss_sim_figDat$amtMiss,1)


# remove "sigma" values, for now
gauss_sim_figDat <- gauss_sim_figDat %>% 
  filter(param != "sigma")

## get the mean and standard deviation for each proportion missing bin
gauss_sim_Means <- gauss_sim_figDat %>%
  group_by(missingness, type, param, round(amtMiss, digits = 1)) %>% 
  summarize(mean_paramDiff = mean(paramDiff_stnd), 
            sd_paramDiff = sd(paramDiff_stnd)) %>% 
  rename("meanAmtMiss" = 'round(amtMiss, digits = 1)')

# make figure 
ggplot(data = gauss_sim_figDat[gauss_sim_figDat$simName %in% c(1:100), ], aes()) + 
  facet_grid(~factor(param, levels = c("intercept", "phi", "light", "discharge")) 
             ~ factor(missingness, levels = c("MAR_lowAutoCor", "MAR_highAutoCor", "MNAR")), scales="free_y") + 
  geom_hline(aes(yintercept = 0), colour = "grey") + 
  geom_point(aes(x = amtMiss, y = paramDiff_stnd, color = as.factor(type)), alpha = .3) +
  #geom_violin(aes(x = factor(amtMiss_bin), y = paramDiff_stnd, color = factor(type)))
  theme_classic() +
  xlab("Proportion of missing data")+ 
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  ylab("standardized parameter estimate [(model param - sim param)/sim param]")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))

ggplot(data = gauss_sim_Means, aes(x = meanAmtMiss, y = mean_paramDiff)) +
  facet_grid(~factor(param, levels = c("intercept", "phi", "light", "discharge")) 
                                                               ~ factor(missingness, levels = c("MAR_lowAutoCor", "MAR_highAutoCor", "MNAR")), scales="free_y") + 
  geom_hline(aes(yintercept = 0), colour = "grey") + 
  geom_errorbar(aes(ymin=mean_paramDiff - sd_paramDiff, ymax=mean_paramDiff + sd_paramDiff, color = as.factor(type)), 
                size=0.3, width=0, position = position_dodge(width=0.03))+
  #geom_ribbon(aes(ymin = mean_paramDiff - sd_paramDiff, ymax = mean_paramDiff + sd_paramDiff, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
  geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
  geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
  theme_classic() +
  xlab("Proportion of missing data")+ 
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  ylab("standardized parameter estimate [(model param - sim param)/sim param]")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))





allmethodsGPP<-ggplot(data=paramallARIMASTAN3, aes(x=as.numeric(missingprop), y=value, color=type))+
  facet_grid(~factor(param, levels=c("intercept", "phi", "light", "discharge"),exclude = NA)~ missingness, scales="free_y")+
  geom_hline(data=trueestdf2, aes(yintercept=value), colour="gray")+
  #geom_hline(data=arimaestdf, aes(yintercept=value), colour="light blue")+
  geom_point(size=0.75, position = position_dodge(width=0.03))+
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), size=0.3, width=0, position = position_dodge(width=0.03))+
  theme_bw()+
  xlab("Proportion of missing data")+ theme(legend.position="top")+theme(legend.title=element_blank())+
  ylab("Parameter estimate")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))





####

