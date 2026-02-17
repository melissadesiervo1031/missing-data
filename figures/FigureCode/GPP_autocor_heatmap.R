## draft of the GPP heat map figure ##

# Load packages
library(tidyverse)
library(ggpubr)

## read in data 

  
gauss_sim_figDat <- readRDS("./data/model_results/gauss_sim_ModelResults.rds")


# remove data for simluation 376... has one really small parameter, which is causing a lot of outliers
gauss_sim_figDat <- gauss_sim_figDat[gauss_sim_figDat$simName != 376,]

##make heat map! ##
# Make heatmaps for Gaussian MCAR data -------------------------------------
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
  #filter(n  > 100) %>% # drop combinations that have fewer than 100 observations
  mutate(tooBigSD = ifelse(paramDiff_SD > 1.96, yes = 1, no = NA)) %>% 
  filter(amtMiss <=.65)
# only consider missingness of 50% or less

# update names for missing data approach
type_names <- c(
  "dropNA_complete" = "Data Deletion-Complete",
  "dropNA_simple" = "Data Deletion-Simple",
  "Kalman Filter" =  "Kalman filter",
  "Multiple Imputations" = "Multiple Imputation",
  "brms" = "Data Augmentation", 
  "Phi" = "Phi", 
  "Beta covariates" = "Beta covariates"
)


figDat2 <- figDat %>% 
  mutate(param = replace(param, param == "Phi", "phi"), 
         param = replace(param, param == "Beta covariates", 'beta'))

## make heatmap for median of parameter recovery
(heatMap_median_MCAR <-ggplot(data = figDat2, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiff_med), size=5) + 
    #scale_fill_viridis_c(name = "value" , option = "A") +
    scale_fill_gradient2(high = "darkblue",
                         mid = "orange",
                         low = "yellow" ,
                         midpoint = -.5, limits = c(-4.8, 2.1),  na.value = "lightgrey") +
    #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    ggh4x::facet_grid2(factor(param, levels = c("phi", "beta")
                               )
                       ~factor(type, levels = c("dropNA_complete", "dropNA_simple", "Kalman Filter", "Multiple Imputations", "brms"), 
                        labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
                      , labeller = label_parsed) +
    xlab("Proportion of missing data")+
    ylab("Autocorrelation in missingness") +
    guides(fill = guide_colorbar("Median Error")) +
    theme_classic() +
    ylim(0,1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Median error of parameter recovery")) 

## make heatmap for SE of parameter recovery
(heatMap_SE_MCAR <-ggplot(data = figDat2, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiffAbsDiff_med), size=5) + 
    #scale_fill_viridis_c(name = "value" , option = "D") +
    scale_fill_gradient2(low= "darkblue",
                         mid = "turquoise",
                         high = "green" ,
                         midpoint = .5, limits = c(0, 0.82),  na.value = "lightgrey") +
    #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    ggh4x::facet_grid2(factor(param, levels = c("phi", "beta")
    )
    ~factor(type, levels = c("dropNA_complete", "dropNA_simple", "Kalman Filter", "Multiple Imputations", "brms"), 
            labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
    , labeller = label_parsed) +
    xlab("Proportion of missing data")+
    ylab("Autocorrelation in missingness") +
    guides(fill = guide_colorbar("Median \nAbsolute Error")) +
    theme_classic() +
    ylim(0,1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Median absolute error of parameter recovery")) 

## make heatmap for SD of parameter recovery
# calculate coverage 
figDat_cov_temp <- gauss_sim_figDat  %>% 
  mutate(CI95_lower = param_value - 1.96*param_se,
         CI95_upper = param_value + 1.96*param_se) %>% 
  filter(param != "sigma") %>% 
  filter(!is.na(param_se))# randomly there are some model runs that don't have SE? 

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
         param = replace(param, param %in% c("discharge", "light"), "beta"),
         #param = replace(param, param == "intercept", "Intercept"),
         param = replace(param, param == "phi", "phi")
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
 
(heatMap_cov_MCAR <- ggplot(data = figDat_covMAR, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=coveragePerc*100), size=5) + 
    #geom_tile(data = figDat_covMAR[figDat_covMAR$coveragePerc >.96,], aes(), fill = "grey") + 
    ggh4x::facet_grid2(factor(param, levels = c("phi", "beta")
    )
    ~factor(type, levels = c("dropNA_complete", "dropNA_simple", "Kalman Filter", "Multiple Imputations", "brms"), 
            labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
    , labeller = label_parsed) +
     #scale_fill_distiller(palette = "Greys", direction = 1, name = "value") +#value = SD
    #scale_fill_viridis_c(begin=0, end=1, option = "G", name = "value" )+ 
    scale_fill_gradient(low = "darkblue",
                        high = "lightblue" ,
                        limits = c(0,100),  na.value = "lightgrey") +
    guides(fill = guide_colorbar("Median \n Absolute Error")) +
    xlab("Proportion of missing data")+
    ylab("Autocorrelation in missingness") +
    guides(fill = guide_colorbar("% Coverage")) +
    theme_classic() +
    ylim(0,1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("% of model runs where the 95% CI includes the simulation parameter"))

## save figures
png(file = "./figures/heatmap_GaussianMCAR_median.png", width = 8, height = 6, units = "in", res = 700)
heatMap_median_MCAR
dev.off()

png(file = "./figures/heatmap_GaussianMCAR_SE.png", width = 8, height = 6, units = "in", res = 700)
heatMap_SE_MCAR
dev.off()

png(file = "./figures/heatmap_GaussianMCAR_coverage.png", width = 8, height = 6, units = "in", res = 700)
heatMap_cov_MCAR
dev.off()

png(file = "./figures/heatmap_GaussianMCAR_all.png", width = 9.5, height = 12, units = "in", res = 700)
ggarrange(heatMap_median_MCAR, heatMap_SE_MCAR, heatMap_cov_MCAR, ncol = 1) %>% 
  annotate_figure(top = text_grob("Parameter recovery from real-valued time series with MCAR data", just = .65, 
                                  size = 16, face = "bold"))
dev.off()



# Get figure of MCAR median error for phi ----------------------------------
figDat3 <- figDat2 %>% 
  filter(param == "phi")
(heatMap_median_MCAR_slice <-ggplot(data = figDat3, aes(x=amtMiss, y=autoCor)) + 
   geom_tile(aes(fill=paramDiff_med), size=5) + 
   #scale_fill_viridis_c(name = "value" , option = "A") +
    scale_fill_gradient2(high = "darkblue",
                         mid = "orange",
                         low = "yellow" ,
                         midpoint = -.5, limits = c(-0.82, .050),  na.value = "lightgrey") +
   #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
   ggh4x::facet_grid2(
   ~factor(type, levels = c("Data Deletion Simple", "Data Deletion CC",  "Multiple imputations","Kalman filter",  "brms"), 
           labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple imputation"', '"Data Augmentation"'))
   , labeller = label_parsed) +
   xlab("Proportion of missing data")+
   ylab("Autocorrelation in missingness") +
   guides(fill = guide_colorbar("Median Error")) +
   theme_classic() +
    ylim(0,1) +
   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
   ggtitle("Missing completely at random: Median error of parameter recovery, \u03D5" )) 


png(file = "./figures/heatmap_GaussianMCAR_justPhi.png", width = 9, height = 3, units = "in", res = 700)
heatMap_median_MCAR_slice
dev.off()

# Make heatmaps for Gaussian MNAR data ------------------------------------
# bin amt missing and autocorr (average paramDiff)
figDat <- gauss_sim_figDat %>% 
  filter(param != "sigma" & 
           param != "intercept" & 
           missingness == "MNAR") %>% 
  mutate(autoCor = 0, 
         amtMiss = round(amtMiss, 1),
         param = replace(param, param %in% c("discharge", "light"), "beta"), # group light and discharge
         param = replace(param, param == "intercept", "Intercept"),
         param = replace(param, param == "phi", "phi")) %>% 
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
    ggh4x::facet_grid2(factor(param, levels = c("phi", "beta")
    )
    ~factor(type, levels = c("Data Deletion CC", "Data Deletion Simple", "Kalman filter", "Multiple imputations", "brms"), 
            labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
    , labeller = label_parsed) +
    scale_fill_viridis_c(begin=1, end=0, option = "A")+
    xlab("Proportion of missing data")+
    ylab("Autocorrelation in missingness") +
    guides(fill = guide_colorbar("Median Error")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Median of parameter error across sims.")) 

## make heatmap for mean of parameter recovery
(heatMap_SE_MNAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiffAbsDiff_med), size=5) + 
    #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    ggh4x::facet_grid2(factor(param, levels = c("phi", "beta")
    )
    ~factor(type, levels = c("Data Deletion CC", "Data Deletion Simple", "Kalman filter", "Multiple imputations", "brms"), 
            labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
    , labeller = label_parsed) +
    scale_fill_viridis_c(begin=1, end=0, option = "D")+
    xlab("Proportion of missing data")+
    ylab("Autocorrelation in missingness") +
    theme_classic() +
    guides(fill = guide_colorbar("Median \n Absolute Error")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Absolute standard error of parameter recovery")) 


figDat_covMNAR$autoCor <- 0
(heatMap_cov_MNAR <- ggplot(data = figDat_covMNAR, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=coveragePerc*100), size=5) + 
    ggh4x::facet_grid2(factor(param, levels = c("phi", "beta")
    )
    ~factor(type, levels = c("Data Deletion CC", "Data Deletion Simple", "Kalman filter", "Multiple imputations", "brms"), 
            labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
    , labeller = label_parsed) +
    #scale_fill_distiller(palette = "Greys", direction = 1, name = "value") +#value = SD
    #scale_fill_viridis_c(begin=0, end=1, option = "G", name = "value" )+ 
    scale_fill_gradient(low = "darkblue",
                        high = "lightblue" ,
                        limits = c(0,100),  na.value = "lightgrey") +
    xlab("Proportion of missing data")+
    guides(fill = guide_colorbar("% Coverage")) +
    ylab("Autocorrelation in missingness") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("% of model runs where the 95% CI \n includes the simulation parameter"))


## save figures
png(file = "./figures/heatmap_GaussianMNAR_median.png", width = 8, height = 6, units = "in", res = 700)
heatMap_medians_MNAR
dev.off()

png(file = "./figures/heatmap_GaussianMNAR_SE.png", width = 8, height = 6, units = "in", res = 700)
heatMap_SE_MNAR
dev.off()

png(file = "./figures/heatmap_GaussianMNAR_coverage.png", width = 8, height = 6, units = "in", res = 700)
heatMap_cov_MNAR
dev.off()

png(file = "./figures/heatmap_GaussianMNAR_all.png", width = 9, height = 12, units = "in", res = 700)
ggarrange(heatMap_medians_MNAR, heatMap_SE_MNAR, heatMap_cov_MNAR, ncol = 1)
dev.off()


