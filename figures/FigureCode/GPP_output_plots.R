#/////////////////
# This script makes figures for GPP simulated data 
# 30 October 2023
#/////////////////

# Load packages
library(tidyverse)
library(ggpubr)
library(forecast)


## read in data 
figDat_temp <- readRDS("./data/model_results/gauss_sim_ModelResults_normPrior.rds")

# remove data for simluation 376... has one really small parameter, which is causing a lot of outliers
figDat_temp <- figDat_temp[figDat_temp$simName != 376,]

# filter for low and high autocor
figDat_temp[figDat_temp$autoCor <=0.3 & !is.na(figDat_temp$autoCor), "missingness"] <- "MAR_lowAutoCor"
figDat_temp[figDat_temp$autoCor > 0.3 & figDat_temp$autoCor < 0.6 &!is.na(figDat_temp$autoCor), "missingness"] <- "MAR_medAutoCor"
figDat_temp[figDat_temp$autoCor  >= 0.6 & !is.na(figDat_temp$autoCor), "missingness"] <- "MAR_highAutoCor"

figDat_lines <- figDat_temp %>% 
  filter(param != "sigma") %>%
  filter(missingness %in% c("MAR_highAutoCor", "MAR_lowAutoCor", "MAR_medAutoCor", "MNAR")) %>% 
  mutate(autoCor = round(autoCor, 1), 
         amtMiss = round(amtMiss, 1)) %>% 
  group_by(missingness, type, param, amtMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) %>% 
  filter(n  > 100)  %>% # drop combinations that have fewer than 300 observations
  filter(amtMiss <=.4)

# Figure of parameter recovery (mean and sd in separate panels) -----------
# figure of means for each model type and level of missingness (with shortened x-axis)
(gauss_sim_MedsFig_trimmed <- ggplot(data = figDat_lines, aes(x = amtMiss, y = paramDiff_med)) +
  facet_grid(~factor(param, levels = c( "intercept","phi", "light", "discharge")) 
                                       ~ factor(missingness, levels = c("MAR_lowAutoCor", "MAR_medAutoCor", "MAR_highAutoCor", "MNAR")),
             scales = "free_y") + 
  geom_hline(aes(yintercept = 0), colour = "grey") + 
  #geom_errorbar(aes(ymin=paramDiff_mean - paramDiff_SD, ymax=paramDiff_mean + paramDiff_SD, color = as.factor(type)), 
                #size=0.3, width=0, position = position_dodge(width=0.03))+
  #geom_ribbon(aes(ymin = paramDiff_mean - paramDiff_SD, ymax = paramDiff_mean + paramDiff_SD, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
  geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
  geom_point(data = data.frame("x" = c(0,.5), "y" = c(-.25, .1)), aes(x = x, y = y), alpha = 0.0000001) + ## add invisible points to change scales
  geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
  theme_classic() +
  xlab("Proportion of missing data")+ 
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  ylab("Median standardized parameter estimate")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
   xlim(c(-0.03,0.43))
  #ylim(c(-0.85,0.55))
)

## make a sub-dataframe that includes points for SDs that are >5 
largeSD <- figDat_lines[figDat_lines$amtMiss <= 0.5 & 
                             figDat_lines$paramDiff_SD > 5,]
  
# figure of SDfor each model type and level of missingness
(gauss_sim_SDFig_trimmed <- ggplot(data = figDat_lines, aes(x = amtMiss, y = paramDiff_SD)) +
  facet_grid(~factor(param, levels = c("intercept","phi", "light", "discharge")) 
             ~ factor(missingness, levels = c("MAR_lowAutoCor", "MAR_medAutoCor", "MAR_highAutoCor", "MNAR"))) + 
  geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
  geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
  geom_hline(aes(yintercept = 0), colour = "grey") + theme_classic() +
  xlab("Proportion of missing data")+ 
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  ylab("SD of standardized parameter estimate")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))+
    xlim(c(-0.03,0.43)) #+ 
    #ylim(c(0,5)) #+
  #geom_point(data = largeSD, aes(x = amtMiss, y = c(3.99,3.99), color = as.factor(type)), 
            # position = position_dodge(width=0.03), pch = 8)
  ) 

# put into one figure
Gauss_paramRecov_trimmed <- ggarrange(gauss_sim_MedsFig_trimmed, gauss_sim_SDFig_trimmed, common.legend = TRUE)

## save results
png(file = "./figures/parameterRecovery_sim_Guassian_meansSD_trimmed.png", width = 9, height = 4, units = "in", res = 700)
Gauss_paramRecov_trimmed
dev.off()


#
figDat_long <- figDat_temp %>% 
  filter(param != "sigma") %>%
  filter(missingness %in% c("MAR_highAutoCor", "MAR_medAutoCor", "MAR_lowAutoCor", "MNAR")) %>% 
  mutate(autoCor = round(autoCor, 1), 
         amtMiss = round(amtMiss, 1)) %>% 
  group_by(missingness, type, param, amtMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) %>% 
  filter(n  > 100)  

# make a figure like the one above, but without a trimmed x axis
# figure of means for each model type and level of missingness
(gauss_sim_MedsFig_reg<- ggplot(data = figDat_long, aes(x = amtMiss, y = paramDiff_med)) +
    facet_grid(~factor(param, levels = c("intercept", "phi", "light", "discharge")) 
               ~ factor(missingness, levels = c("MAR_lowAutoCor", "MAR_medAutoCor", "MAR_highAutoCor", "MNAR"))) + 
    geom_hline(aes(yintercept = 0), colour = "grey") + 
    #geom_errorbar(aes(ymin=paramDiff_mean - paramDiff_SD, ymax=paramDiff_mean + paramDiff_SD, color = as.factor(type)), 
    #size=0.3, width=0, position = position_dodge(width=0.03))+
    #geom_ribbon(aes(ymin = paramDiff_mean - paramDiff_SD, ymax = paramDiff_mean + paramDiff_SD, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("Mean standardized parameter estimate")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) 
)

# figure of SD for each model type and level of missingness
(gauss_sim_SDFig_reg<- ggplot(data = figDat_long, aes(x = amtMiss, y = paramDiff_SD)) +
    facet_grid(~factor(param, levels = c("intercept", "phi", "light", "discharge")) 
               ~ factor(missingness, levels = c("MAR_lowAutoCor", "MAR_medAutoCor", "MAR_highAutoCor", "MNAR"))) + 
    geom_hline(aes(yintercept = 0), colour = "grey") + 
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    #geom_ribbon(aes(ymin = 0, ymax = paramDiff_SD, fill = as.factor(type), color = as.factor(type)), alpha = .2, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("SD of standardized parameter estimate")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))
) 
# put into one figure
Gauss_paramRecov <- ggarrange(gauss_sim_MedsFig_reg, gauss_sim_SDFig_reg, common.legend = TRUE)

## save results
png(file = "./figures/parameterRecovery_sim_Guassian_medsSD.png", width = 9, height = 4, units = "in", res = 700)
Gauss_paramRecov
dev.off()
# 
# 
# # 95% CI Error bar plots to show spread of complete parameter recovery data  --------
# (ErrorBarPlots <- ggplot(data = figDat_long, aes(x = amtMiss, y = paramDiff_med)) +
#    facet_grid(~factor(param, levels = c("phi", "intercept",  "light", "discharge")) 
#               ~ factor(missingness, levels = c("MAR_lowAutoCor", "MAR_highAutoCor", "MNAR"))) + 
#    geom_hline(aes(yintercept = 0), colour = "grey") + 
#    geom_errorbar(aes(ymin=paramDiff_med - 1.96*paramDiff_SD, ymax=paramDiff_med + 1.96*paramDiff_SD, color = as.factor(type)), 
#       size=0.3, width=0, position = position_dodge(width=0.07))+
#    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.07)) +
#    theme_classic() +
#    xlab("Proportion of missing data")+ 
#    theme(legend.position="top")+
#    theme(legend.title=element_blank())+
#    ylab("Mean standardized parameter estimate")+ 
#    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))
# )
# 
# 
# figDat_all <- figDat_temp %>% 
#   filter(param != "sigma") %>%
#   filter(missingness %in% c("MAR_highAutoCor", "MAR_lowAutoCor", "MNAR")) %>% 
#   mutate(autoCor = round(autoCor, 1), 
#          amtMiss = round(amtMiss, 1)) %>% 
#   filter(amtMiss <= 0.5)
# 
# (gauss_sim_violin <- ggplot(data = figDat_all) +
#   facet_grid(~factor(param, levels = c("intercept", "phi", "light", "discharge")) 
#              ~ factor(missingness, levels = c("MAR_lowAutoCor", "MAR_,medAutoCor", "MAR_highAutoCor", "MNAR"))) + 
#   geom_hline(aes(yintercept = 0), colour = "grey") + 
#   geom_violin(aes(x = as.factor(amtMiss), y = paramDiff, color = type)) +
#   theme_classic() +
#   xlab("Proportion of missing data")+ 
#   theme(legend.position="top")+
#   theme(legend.title=element_blank())+
#   ylab("Mean standardized parameter estimate")+ 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)))
# 
# 
# ## save figure
# png(file = "./figures/parameterRecovery_sim_Guassian_95CIs.png", width = 9, height = 4, units = "in", res = 700)
# ErrorBarPlots
# dev.off()


## trying to figure out brms phi issue
phi <- figDat_temp[figDat_temp$param == "phi",]
# calculate the unstandardized difference of values
phi <- phi %>% 
  mutate(paramDiff_actual = value - param_simVal)

phi_means <- phi %>% 
  filter(param != "sigma") %>%
  filter(missingness %in% c("MAR_highAutoCor", "MAR_lowAutoCor", "MNAR")) %>% 
  mutate(autoCor = round(autoCor, 1), 
         amtMiss = round(amtMiss, 1)) %>% 
  group_by(missingness, type, param, amtMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff),
            paramDiff_actual_mean = mean(paramDiff_actual, na.rm = TRUE),
            param_actual_sd = sd(paramDiff_actual, na.rm = TRUE)) %>% 
  filter(n  > 100)  %>% # drop combinations that have fewer than 300 observations
  filter(amtMiss <=.5)
  
# calculate quantiles
phi_quants <- phi %>% 
  group_by(type) %>% 
  summarize("quantile_value" = quantile(x = paramDiff, probs = c(0.01, .25, .5, .75, .99)),
            "quantile_level" = c(.01, .25, .5, .75, .99))

  
ggplot(data = phi) +
  geom_histogram(aes(paramDiff, col = type)) + 
  facet_wrap(~type)

ggplot(data = phi) + 
  geom_violin(aes(x = type, y = paramDiff, col = type)) +
  xlab("Missingness Approach") + 
  ylab("[(true param-sim param)/abs(sim param)]") +
  ggtitle("phi parameter recovery")

ggplot(data = phi) + 
  geom_vline(aes(xintercept = c(0)), col = "grey", lty = 2) +
  geom_vline(aes(xintercept = c(1)), col = "grey", lty = 2) +
  geom_histogram(aes(param_simVal), col = "grey20", alpha = .3) +
  geom_histogram(aes(value, col = type, fill = type), alpha = .5) + 
  facet_wrap(~type) + 
  xlab("model estimate of phi") + 
  theme_classic()
  
ggplot(data = phi) + 
  geom_vline(aes(xintercept = c(0)), col = "grey", lty = 2) +
  geom_vline(data = phi_quants[phi_quants$quantile_level == .01,], aes(xintercept =quantile_value)) + 
  geom_vline(data = phi_quants[phi_quants$quantile_level == .99,], aes(xintercept =quantile_value)) + 
  geom_histogram(aes(paramDiff, col = type, fill = type), alpha = .5) + 
  facet_wrap(~type) + 
  xlab("standardized parameter differences ") + 
  theme_classic() + 
  xlim(c(min(phi_quants$quantile_value), max(phi_quants$quantile_value)))


