## draft of the GPP heat map figure ##

# Load packages
library(tidyverse)
library(ggpubr)

## read in data 
gauss_sim_ModelResults <- readRDS("./data/model_results/gauss_sim_ModelResults.rds")
gauss_sim_ModelResults <- unique(gauss_sim_ModelResults)

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
# 
gauss_sim_ModelResults[gauss_sim_ModelResults$missingness == "MAR" & 
                                 is.na(gauss_sim_ModelResults$autoCor), "autoCor"] <- 0

# fix values for MNAR (remove autocor values)
gauss_sim_ModelResults[gauss_sim_ModelResults$missingness == "MNAR", "autoCor"] <- NA

# remove values for models fitted to time series with no missingness (Doesn't work for all model approaches)
gauss_sim_figDat <- gauss_sim_ModelResults[gauss_sim_ModelResults$missingprop_autocor != "y_noMiss",]

simDat <-
  gauss_sim_figDat %>%
  select(simName, phi_sim, intercept_sim, light_sim, discharge_sim) %>%
  pivot_longer(cols = c(phi_sim, intercept_sim, light_sim, discharge_sim),
               names_to = "param",
               values_to = "param_simVal",
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>%
  unique()


# remove columns for simulation data
gauss_sim_figDat <- gauss_sim_figDat %>% 
  select(-phi_sim, -intercept_sim, -light_sim, -discharge_sim)

gauss_sim_figDat <- gauss_sim_figDat %>% 
  left_join(simDat, by = c("simName", "param"))

# calculate the standardized difference between parameter estimates and simulated values??
gauss_sim_figDat <- gauss_sim_figDat %>% 
  dplyr::mutate(paramDiff = ((value - param_simVal)/abs(param_simVal)))

##make heat map! ##
# Make heatmaps for Gaussian MAR data -------------------------------------
# bin amt missing and autocorr (average paramDiff)
figDat <- gauss_sim_figDat %>% 
  filter(param != "sigma" & 
           missingness == "MAR") %>% 
  mutate(autoCor = round(autoCor, 1), 
         amtMiss = round(amtMiss, 1)) %>% 
  group_by(missingness, type, param, autoCor, amtMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) %>% 
  filter(n  > 300) %>% # drop combinations that have fewer than 300 observations
  mutate(tooBigSD = ifelse(paramDiff_SD > 1.96, yes = 1, no = NA)) 

## make heatmap for mean of parameter recovery
(heatMap_mean_MAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiff_mean), size=5) + 
    scale_fill_viridis_c(name = "value" ) +
    facet_grid(~factor(figDat$param, levels = c("phi", "intercept","light", "discharge")) ~ type) +
    #scale_fill_viridis_c(begin=1, end=0, option = "inferno")+
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Mean of standardized parameter estimates")) 

## make heatmap for SD of parameter recovery
(heatMap_SD_MAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiff_SD), size=5) + 
    facet_grid(~factor(figDat$param, levels = c("phi", "intercept","light", "discharge")) ~ type) +
    scale_fill_viridis_c(begin=1, end=0, option = "plasma", name = "value" )+
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("SD of standardized parameter estimates"))

## save figures
png(file = "./figures/heatmap_GaussianMAR_mean.png", width = 7, height = 6, units = "in", res = 700)
heatMap_mean_MAR
dev.off()

png(file = "./figures/heatmap_GaussianMAR_SD.png", width = 7, height = 6, units = "in", res = 700)
heatMap_SD_MAR
dev.off()

png(file = "./figures/heatmap_GaussianMAR_all.png", width = 12.5, height = 6, units = "in", res = 700)
ggarrange(heatMap_mean_MAR, heatMap_SD_MAR)
dev.off()

# Make heatmaps for Gaussian MNAR data ------------------------------------
# bin amt missing and autocorr (average paramDiff)
figDat <- gauss_sim_figDat %>% 
  filter(param != "sigma" & 
           missingness == "MNAR") %>% 
  mutate(autoCor = 0, 
         amtMiss = round(amtMiss, 1)) %>% 
  group_by(missingness, type, param, autoCor, amtMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) %>% 
  filter(n  > 300) %>% # drop combinations that have fewer than 300 observations
  mutate(tooBigSD = ifelse(paramDiff_SD > 1.96, yes = 1, no = NA)) 

## make heatmap for mean of parameter recovery
(heatMap_mean_MNAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiff_mean), size=5) + 
    scale_fill_viridis_c(name = "value") +
    facet_grid(~factor(figDat$param, levels = c("phi", "intercept","light", "discharge")) ~ type) +
    #scale_fill_viridis_c(begin=1, end=0, option = "inferno")+
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Mean of standardized parameter estimates")) 

## make heatmap for SD of parameter recovery
(heatMap_SD_MNAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiff_SD), size=5) + 
    facet_grid(~factor(figDat$param, levels = c("phi", "intercept","light", "discharge")) ~ type) +
    scale_fill_viridis_c(begin=1, end=0, option = "plasma", name = "value" )+
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("SD of standardized parameter estimates"))

## save figures
png(file = "./figures/heatmap_GaussianMNAR_mean.png", width = 7, height = 6, units = "in", res = 700)
heatMap_mean_MNAR
dev.off()

png(file = "./figures/heatmap_GaussianMNAR_SD.png", width = 7, height = 6, units = "in", res = 700)
heatMap_SD_MNAR
dev.off()

png(file = "./figures/heatmap_GaussianMNAR_all.png", width = 12.5, height = 6, units = "in", res = 700)
ggarrange(heatMap_mean_MNAR, heatMap_SD_MNAR)
dev.off()
