#/////////////////
# This script makes figures for GPP simulated data 
# 30 October 2023
#/////////////////

# Load packages
library(here)
library(tidyverse)
library(ggpubr)

# read in data resulting from beartooth
gauss_sim_ModelResults <- readRDS("./data/model_results/gauss_sim_ModelResults.rds")

# make columns for "autocor" and "missingness"
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
         | (missingness == "MAR" & autoCor <= 0.2) # get low autocor MAR data
           | (missingness == "MAR" & autoCor >= 0.8 )) 

gauss_sim_figDat[gauss_sim_figDat$autoCor <=0.2 & !is.na(gauss_sim_figDat$autoCor), "missingness"] <- "MAR_lowAutoCor"
gauss_sim_figDat[gauss_sim_figDat$autoCor  >= 0.8 & !is.na(gauss_sim_figDat$autoCor), "missingness"] <- "MAR_highAutoCor"

# remove values for models fitted to time series with no missingness (Doesn't work for all model approaches)
gauss_sim_figDat <- gauss_sim_figDat[gauss_sim_figDat$missingprop_autocor != "y_noMiss",]

# calculate the standardized difference between parameter estimates and simulated values??
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

# remove duplicate rows that got introduced to BRMS mnar portion somehow
gauss_sim_figDat <- unique(gauss_sim_figDat)

## get the mean and standard deviation for each proportion missing bin
gauss_sim_Means <- gauss_sim_figDat %>%
  group_by(missingness, type, param, round(amtMiss, digits = 1)) %>% 
  summarize(mean_paramDiff = mean(paramDiff_stnd), 
            sd_paramDiff = sd(paramDiff_stnd),
            n_bin = length(paramDiff_stnd)) %>% 
  rename("meanAmtMiss" = 'round(amtMiss, digits = 1)')


# Figure of parameter recovery (mean and sd in separate panels) -----------
# figure of means for each model type and level of missingness (with shortened x-axis)
(gauss_sim_MeansFig_trimmed <- ggplot(data = gauss_sim_Means, aes(x = meanAmtMiss, y = mean_paramDiff)) +
  facet_grid(~factor(param, levels = c("intercept", "phi", "light", "discharge")) 
                                       ~ factor(missingness, levels = c("MAR_lowAutoCor", "MAR_highAutoCor", "MNAR"))) + 
  geom_hline(aes(yintercept = 0), colour = "grey") + 
  #geom_errorbar(aes(ymin=mean_paramDiff - sd_paramDiff, ymax=mean_paramDiff + sd_paramDiff, color = as.factor(type)), 
                #size=0.3, width=0, position = position_dodge(width=0.03))+
  #geom_ribbon(aes(ymin = mean_paramDiff - sd_paramDiff, ymax = mean_paramDiff + sd_paramDiff, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
  geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
  geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
  theme_classic() +
  xlab("Proportion of missing data")+ 
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  ylab("Mean standardized parameter estimate")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
   xlim(c(-0.05,.55)) +
   ylim(c(-0.85,0.8))
)

## make a sub-dataframe that includes points for SDs that are >5 
largeSD <- gauss_sim_Means[gauss_sim_Means$meanAmtMiss <= 0.5 & 
                             gauss_sim_Means$sd_paramDiff > 4,]
  
# figure of means for each model type and level of missingness
(gauss_sim_SDFig_trimmed <- ggplot(data = gauss_sim_Means, aes(x = meanAmtMiss, y = sd_paramDiff)) +
  facet_grid(~factor(param, levels = c("intercept", "phi", "light", "discharge")) 
             ~ factor(missingness, levels = c("MAR_lowAutoCor", "MAR_highAutoCor", "MNAR"))) + 
  geom_hline(aes(yintercept = 1.96), colour = "grey", lty = 2) + 
  geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
  geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
  #geom_ribbon(aes(ymin = 0, ymax = sd_paramDiff, fill = as.factor(type), color = as.factor(type)), alpha = .2, position = position_dodge(width=0.03)) +
  theme_classic() +
  xlab("Proportion of missing data")+ 
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  ylab("SD of standardized parameter estimate")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))+
    xlim(c(-0.05,.55)) + 
    ylim(c(0,4)) +
  geom_point(data = largeSD, aes(x = meanAmtMiss, y = c(3.99,3.99,3.99,3.99,3.99,3.99), color = as.factor(type)), 
             position = position_dodge(width=0.03), pch = 8)
  ) 

# put into one figure
Gauss_paramRecov_trimmed <- ggarrange(gauss_sim_MeansFig_trimmed, gauss_sim_SDFig_trimmed, common.legend = TRUE)

## save results
png(file = "./figures/parameterRecovery_sim_Guassian_meansSD_trimmed.png", width = 9, height = 4, units = "in", res = 700)
Gauss_paramRecov_trimmed
dev.off()


## make a figure like the one above, but without a trimmed x axis
# figure of means for each model type and level of missingness
(gauss_sim_MeansFig_reg<- ggplot(data = gauss_sim_Means, aes(x = meanAmtMiss, y = mean_paramDiff)) +
    facet_grid(~factor(param, levels = c("intercept", "phi", "light", "discharge")) 
               ~ factor(missingness, levels = c("MAR_lowAutoCor", "MAR_highAutoCor", "MNAR"))) + 
    geom_hline(aes(yintercept = 0), colour = "grey") + 
    #geom_errorbar(aes(ymin=mean_paramDiff - sd_paramDiff, ymax=mean_paramDiff + sd_paramDiff, color = as.factor(type)), 
    #size=0.3, width=0, position = position_dodge(width=0.03))+
    #geom_ribbon(aes(ymin = mean_paramDiff - sd_paramDiff, ymax = mean_paramDiff + sd_paramDiff, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("Mean standardized parameter estimate")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) 
)
# figure of means for each model type and level of missingness
(gauss_sim_SDFig_reg<- ggplot(data = gauss_sim_Means, aes(x = meanAmtMiss, y = sd_paramDiff)) +
    facet_grid(~factor(param, levels = c("intercept", "phi", "light", "discharge")) 
               ~ factor(missingness, levels = c("MAR_lowAutoCor", "MAR_highAutoCor", "MNAR"))) + 
    geom_hline(aes(yintercept = 1.96), colour = "grey", lty = 2) + 
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    #geom_ribbon(aes(ymin = 0, ymax = sd_paramDiff, fill = as.factor(type), color = as.factor(type)), alpha = .2, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("SD of standardized parameter estimate")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))
) 
# put into one figure
Gauss_paramRecov <- ggarrange(gauss_sim_MeansFig_reg, gauss_sim_SDFig_reg, common.legend = TRUE)

## save results
png(file = "./figures/parameterRecovery_sim_Guassian_meansSD.png", width = 9, height = 4, units = "in", res = 700)
Gauss_paramRecov
dev.off()


# 95% CI Error bar plots to show spread of complete parameter recovery data  --------
(ErrorBarPlots <- ggplot(data = gauss_sim_Means, aes(x = meanAmtMiss, y = mean_paramDiff)) +
   facet_grid(~factor(param, levels = c("intercept", "phi", "light", "discharge")) 
              ~ factor(missingness, levels = c("MAR_lowAutoCor", "MAR_highAutoCor", "MNAR")), scales = "free_y") + 
   geom_hline(aes(yintercept = 0), colour = "grey") + 
   geom_errorbar(aes(ymin=mean_paramDiff - 1.96*sd_paramDiff, ymax=mean_paramDiff + 1.96*sd_paramDiff, color = as.factor(type)), 
   size=0.3, width=0, position = position_dodge(width=0.07))+
   geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.07)) +
   theme_classic() +
   xlab("Proportion of missing data")+ 
   theme(legend.position="top")+
   theme(legend.title=element_blank())+
   ylab("Mean standardized parameter estimate")+ 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))
)
## save figure
png(file = "./figures/.png", width = 9, height = 4, units = "in", res = 700)
ErrorBarPlots
dev.off()
