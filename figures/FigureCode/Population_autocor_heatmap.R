## draft of the GPP heat map figure ##

# Load packages
library(tidyverse)
library(ggpubr)

## read in data 
ric_sim_figDat <- readRDS(file = "./data/model_results/ricker_sim_ModelResultsLong.rds")

##make heat map! ##
# Make heatmaps for Gaussian MAR data -------------------------------------
# bin amt missing and autocorr (average paramDiff)
figDat <- ric_sim_figDat %>%
  mutate(autoCorr = round(autoCorr, 1), 
         propMiss = round(propMiss, 1)) %>% 
  group_by(missingness, type, param, autoCorr, propMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) %>% 
  #filter(n  > 100) %>% # drop combinations that have fewer than 100 observations
  filter(propMiss <=.5)
# only consider missingness of 50% or less

## make heatmap for mean of parameter recovery
(heatMap_median_MAR <-ggplot(data = figDat, aes(x=propMiss, y=autoCorr)) + 
    geom_tile(aes(fill=paramDiff_med), size=5) + 
    #scale_fill_viridis_c(name = "value" , option = "magma") +
    scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    facet_grid(~factor(figDat$param, levels = c("alpha", "r")) ~ type) +
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Median of parameter bias across sims.")) 

## make heatmap for SD of parameter recovery
(heatMap_SD_MAR <-ggplot(data = figDat, aes(x=propMiss, y=autoCorr)) + 
    geom_tile(aes(fill=paramDiff_SD), size=5) + 
    facet_grid(~factor(figDat$param, levels = c("alpha", "r")) ~ type) +
    scale_fill_distiller(palette = "Greys", direction = 1, name = "value") +
    #viridis_c(begin=1, end=0, option = "plasma", name = "value" )+
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("SD of parameter bias across sims."))

## save figures
png(file = "./figures/heatmap_PoissonMAR_median.png", width = 7, height = 6, units = "in", res = 700)
heatMap_median_MAR
dev.off()

png(file = "./figures/heatmap_PoissonMAR_SD.png", width = 7, height = 6, units = "in", res = 700)
heatMap_SD_MAR
dev.off()

png(file = "./figures/heatmap_PoissonMAR_all.png", width = 12.5, height = 6, units = "in", res = 700)
ggarrange(heatMap_median_MAR, heatMap_SD_MAR)
dev.off()
