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
           param != "intercept" &
           missingness == "MAR") %>% 
  mutate(autoCor = round(autoCor, 1), 
         amtMiss = round(amtMiss, 1),
         param = replace(param, param %in% c("discharge", "light"), "Beta covariates"), # group light and discharge
         param = replace(param, param == "intercept", "Intercept"),
         param = replace(param, param == "phi", "Phi")) %>% 
  group_by(missingness, type, param, autoCor, amtMiss) %>% 
  summarize(paramDiffAbsDiff_mean = mean(paramDiff_absDiff, na.rm = TRUE),
            paramDiffAbsDiff_med = median(paramDiff_absDiff, na.rm = TRUE),
            paramDiffAbsDiff_SD = sd(paramDiff_absDiff, na.rm = TRUE),
            n_paramDiffAbsDiff = length(paramDiff_absDiff), 
            paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) %>% 
  filter(n  > 100) %>% # drop combinations that have fewer than 100 observations
  mutate(tooBigSD = ifelse(paramDiff_SD > 1.96, yes = 1, no = NA)) %>% 
  filter(amtMiss <=.5)
# only consider missingness of 50% or less

# update names for missing data approach
type_names <- c(
  "Data Deletion CC" = "Data Deletion-Complete",
  "Data Deletion Simple" = "Data Deletion-Simple",
  "Kalman filter" =  "Kalman filter",
  "Multiple imputations" = "Multiple imputations",
  "brms" = "Data Augmentation", 
  "Phi" = "Phi", 
  "Beta covariates" = "Beta covariates"
)

## make heatmap for median of parameter recovery
(heatMap_median_MAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiff_med), size=5) + 
    scale_fill_viridis_c(name = "value" , option = "A") +
    #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    facet_grid(~factor(figDat$param, levels = c("Phi", "Beta covariates")) ~ type, 
               labeller = as_labeller(type_names)) +
    xlab("Proportion of missing data")+
    ylab("Autocorrelation in missingness") +
    guides(fill = guide_colorbar("Bias")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Median of parameter bias across sims.")) 

## make heatmap for SE of parameter recovery
(heatMap_SE_MAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiffAbsDiff_med), size=5) + 
    scale_fill_viridis_c(name = "value" , option = "D") +
    #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    facet_grid(~factor(figDat$param, levels = c("Phi", "Beta covariates")) ~ type, 
               labeller = as_labeller(type_names)) +xlab("Proportion of missing data")+
    ylab("Autocorrelation in missingness") +
    guides(fill = guide_colorbar("Stand. Error")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Absolute Standard Error of Parameter Recovery")) 

## make heatmap for SD of parameter recovery
# calculate coverage 
figDat_cov_temp <- gauss_sim_figDat  %>% 
  mutate(CI95_lower = value - 1.96*SE,
         CI95_upper = value + 1.96*SE) %>% 
  filter(param != "sigma") %>% 
  filter(!is.na(SE))# randomly there are some model runs that don't have SE? 

# is the true parameter within the 95% CI? 
figDat_cov_temp$coverage <- c(figDat_cov_temp$param_simVal >= figDat_cov_temp$CI95_lower & 
                                    figDat_cov_temp$param_simVal <= figDat_cov_temp$CI95_upper)

## count the # of models w/ and without coverage for each bin of missingness and autocorrelation
figDat_cov <- figDat_cov_temp %>% 
  filter(param != "sigma",
         param != "intercept",
         amtMiss <=.5,) %>% 
  mutate(autoCor = round(autoCor, 1), 
         amtMiss = round(amtMiss, 1),
         param = replace(param, param %in% c("discharge", "light"), "Beta covariates"),
         #param = replace(param, param == "intercept", "Intercept"),
         param = replace(param, param == "phi", "Phi")
  ) %>% 
  group_by(missingness, type, param, autoCor, amtMiss) %>% 
  summarize(coverageNumber = sum(coverage), # the number of models that have coverage
            modelRunN = length(!is.na(coverage))# the total number of models 
            ) %>% 
  mutate(coveragePerc = coverageNumber/modelRunN)
    
figDat_covMAR <- figDat_cov %>% 
  filter(missingness == "MAR")
figDat_covMNAR <- figDat_cov %>% 
  filter(missingness %in% c("MNAR"))
 
(heatMap_cov_MAR <- ggplot(data = figDat_covMAR, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=coveragePerc*100), size=5) + 
    facet_grid(~factor(figDat_covMAR$param, levels = c("Phi", "Beta covariates")) ~ type, 
               labeller = as_labeller(type_names)) +
     #scale_fill_distiller(palette = "Greys", direction = 1, name = "value") +#value = SD
    scale_fill_viridis_c(begin=0, end=1, option = "G", name = "value" )+ 
    guides(fill = guide_colorbar("Stand. Error")) +
    xlab("Proportion of missing data")+
    ylab("Autocorrelation in missingness") +
    guides(fill = guide_colorbar("% Coverage")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("% of model runs where the 95% CI \n includes the simulation parameter"))

## save figures
png(file = "./figures/heatmap_GaussianMAR_median.png", width = 7, height = 6, units = "in", res = 700)
heatMap_median_MAR
dev.off()

png(file = "./figures/heatmap_GaussianMAR_SE.png", width = 7, height = 6, units = "in", res = 700)
heatMap_SE_MAR
dev.off()

png(file = "./figures/heatmap_GaussianMAR_coverage.png", width = 7, height = 6, units = "in", res = 700)
heatMap_cov_MAR
dev.off()

png(file = "./figures/heatmap_GaussianMAR_all.png", width = 12.5, height = 6, units = "in", res = 700)
ggarrange(heatMap_median_MAR, heatMap_SE_MAR, heatMap_cov_MAR)
dev.off()

# Make heatmaps for Gaussian MNAR data ------------------------------------
# bin amt missing and autocorr (average paramDiff)
figDat <- gauss_sim_figDat %>% 
  filter(param != "sigma" & 
           param != "intercept" & 
           missingness == "MNAR") %>% 
  mutate(autoCor = 0, 
         amtMiss = round(amtMiss, 1),
         param = replace(param, param %in% c("discharge", "light"), "Beta covariates"), # group light and discharge
         param = replace(param, param == "intercept", "Intercept"),
         param = replace(param, param == "phi", "Phi")) %>% 
  group_by(missingness, type, param, autoCor, amtMiss) %>% 
  summarize(paramDiffAbsDiff_mean = mean(paramDiff_absDiff, na.rm = TRUE),
            paramDiffAbsDiff_med = median(paramDiff_absDiff, na.rm = TRUE),
            paramDiffAbsDiff_SD = sd(paramDiff_absDiff, na.rm = TRUE),
            n_paramDiffAbsDiff = length(paramDiff_absDiff),
    paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) %>% 
  filter(n  > 300) %>% # drop combinations that have fewer than 300 observations
  mutate(tooBigSD = ifelse(paramDiff_SD > 1.96, yes = 1, no = NA)) 

## make heatmap for mean of parameter recovery
(heatMap_medians_MNAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiff_med), size=5) + 
    #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    facet_grid(~factor(figDat$param, levels = c("Phi", "Beta covariates")) ~ type, 
               labeller = as_labeller(type_names)) +
    scale_fill_viridis_c(begin=1, end=0, option = "A")+
    xlab("Proportion of missing data")+
    ylab("Autocorrelation in missingness") +
    guides(fill = guide_colorbar("Bias")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Median of parameter bias across sims.")) 

## make heatmap for mean of parameter recovery
(heatMap_SE_MNAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiffAbsDiff_med), size=5) + 
    #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    facet_grid(~factor(figDat$param, levels = c("Phi", "Beta covariates")) ~ type,
               labeller = as_labeller(type_names)) +
    scale_fill_viridis_c(begin=1, end=0, option = "D")+
    xlab("Proportion of missing data")+
    ylab("Autocorrelation in missingness") +
    theme_classic() +
    guides(fill = guide_colorbar("Stand. Error")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Absolute standard error of parameter recovery")) 


figDat_covMNAR$autoCor <- 0
(heatMap_cov_MNAR <- ggplot(data = figDat_covMNAR, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=coveragePerc*100), size=5) + 
    facet_grid(~factor(figDat_covMNAR$param, levels = c("Phi", "Beta covariates")) ~ type,
               labeller = as_labeller(type_names)) +
    #scale_fill_distiller(palette = "Greys", direction = 1, name = "value") +#value = SD
    scale_fill_viridis_c(begin=0, end=1, option = "G", name = "value" )+ 
    xlab("Proportion of missing data")+
    guides(fill = guide_colorbar("% Coverage")) +
    ylab("Autocorrelation in missingness") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("% of model runs where the 95% CI \n includes the simulation parameter"))


## save figures
png(file = "./figures/heatmap_GaussianMNAR_median.png", width = 7, height = 6, units = "in", res = 700)
heatMap_medians_MNAR
dev.off()

png(file = "./figures/heatmap_GaussianMNAR_SE.png", width = 7, height = 6, units = "in", res = 700)
heatMap_SE_MNAR
dev.off()

png(file = "./figures/heatmap_GaussianMNAR_coverage.png", width = 7, height = 6, units = "in", res = 700)
heatMap_cov_MNAR
dev.off()

png(file = "./figures/heatmap_GaussianMNAR_all.png", width = 12.5, height = 6, units = "in", res = 700)
ggarrange(heatMap_medians_MNAR, heatMap_SE_MNAR, heatMap_cov_MNAR)
dev.off()
