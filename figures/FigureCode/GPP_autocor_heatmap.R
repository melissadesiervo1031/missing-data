## draft of the GPP heat map figure ##

# Load packages
library(tidyverse)
library(ggpubr)

## read in data 
gauss_sim_figDat <- readRDS("./data/model_results/gauss_sim_ModelResults_normPrior.rds")

# remove data for simluation 376... has one really small parameter, which is causing a lot of outliers
gauss_sim_figDat <- gauss_sim_figDat[gauss_sim_figDat$simName != 376,]

##make heat map! ##
# Make heatmaps for Gaussian MAR data -------------------------------------
# bin amt missing and autocorr (average paramDiff)
figDat <- gauss_sim_figDat %>% 
  filter(param != "sigma" & 
           missingness == "MAR") %>% 
  mutate(autoCor = round(autoCor, 1), 
         amtMiss = round(amtMiss, 1),
         param = replace(param, param %in% c("discharge", "light"), "Beta covariates"), # group light and discharge
         param = replace(param, param == "intercept", "Intercept"),
         param = replace(param, param == "phi", "Phi")) %>% 
  group_by(missingness, type, param, autoCor, amtMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) %>% 
  filter(n  > 100) %>% # drop combinations that have fewer than 100 observations
  mutate(tooBigSD = ifelse(paramDiff_SD > 1.96, yes = 1, no = NA)) %>% 
  filter(amtMiss <=.5)
# only consider missingness of 50% or less

## make heatmap for mean of parameter recovery
(heatMap_median_MAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiff_med), size=5) + 
    #scale_fill_viridis_c(name = "value" , option = "magma") +
    scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    facet_grid(~factor(figDat$param, levels = c("Intercept","Phi", "Beta covariates")) ~ type) +
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Median of parameter bias across sims.")) 

## make heatmap for SD of parameter recovery
(heatMap_SD_MAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiff_SD), size=5) + 
    facet_grid(~factor(figDat$param, levels = c("Intercept","Phi", "Beta covariates")) ~ type) +
    scale_fill_distiller(palette = "Greys", direction = 1, name = "value") +
  #viridis_c(begin=1, end=0, option = "plasma", name = "value" )+
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("SD of parameter bias across sims."))

## save figures
png(file = "./figures/heatmap_GaussianMAR_median.png", width = 7, height = 6, units = "in", res = 700)
heatMap_median_MAR
dev.off()

png(file = "./figures/heatmap_GaussianMAR_SD.png", width = 7, height = 6, units = "in", res = 700)
heatMap_SD_MAR
dev.off()

png(file = "./figures/heatmap_GaussianMAR_all.png", width = 12.5, height = 6, units = "in", res = 700)
ggarrange(heatMap_median_MAR, heatMap_SD_MAR)
dev.off()

# Make heatmaps for Gaussian MNAR data ------------------------------------
# bin amt missing and autocorr (average paramDiff)
figDat <- gauss_sim_figDat %>% 
  filter(param != "sigma" & 
           missingness == "MNAR") %>% 
  mutate(autoCor = 0, 
         amtMiss = round(amtMiss, 1),
         param = replace(param, param %in% c("discharge", "light"), "Beta covariates"), # group light and discharge
         param = replace(param, param == "intercept", "Intercept"),
         param = replace(param, param == "phi", "Phi")) %>% 
  group_by(missingness, type, param, autoCor, amtMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) %>% 
  filter(n  > 300) %>% # drop combinations that have fewer than 300 observations
  mutate(tooBigSD = ifelse(paramDiff_SD > 1.96, yes = 1, no = NA)) 

## make heatmap for mean of parameter recovery
(heatMap_medians_MNAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiff_med), size=5) + 
    scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    facet_grid(~factor(figDat$param, levels = c("Intercept","Phi", "Beta covariates")) ~ type) +
    #scale_fill_viridis_c(begin=1, end=0, option = "inferno")+
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Mean of standardized parameter estimates")) 

## make heatmap for SD of parameter recovery
(heatMap_SD_MNAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiff_SD), size=5) + 
    scale_fill_distiller(palette = "Greys", direction = 1, name = "value")+
    facet_grid(~factor(figDat$param, levels = c("Intercept","Phi", "Beta covariates")) ~ type) +
    #scale_fill_viridis_c(begin=1, end=0, option = "plasma", name = "value" )+
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("SD of standardized parameter estimates"))

## save figures
png(file = "./figures/heatmap_GaussianMNAR_median.png", width = 7, height = 6, units = "in", res = 700)
heatMap_medians_MNAR
dev.off()

png(file = "./figures/heatmap_GaussianMNAR_SD.png", width = 7, height = 6, units = "in", res = 700)
heatMap_SD_MNAR
dev.off()

png(file = "./figures/heatmap_GaussianMNAR_all.png", width = 12.5, height = 6, units = "in", res = 700)
ggarrange(heatMap_medians_MNAR, heatMap_SD_MNAR)
dev.off()
