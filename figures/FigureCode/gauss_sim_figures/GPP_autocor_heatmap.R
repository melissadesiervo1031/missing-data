## draft of the GPP heat map figure ##

# Load packages
library(tidyverse)
library(ggpubr)

## read in data 
## read in data 
figDat_temp <- readRDS("./data/model_results/gauss_sim_ModelResults.rds") ## read in file from dropbox..too big for Githu##


# remove data for simluation 376... has one really small parameter, which is causing a lot of outliers

#figDat_temp <- figDat_temp[!(figDat_temp$simName %in% c(376, 831, 816, 461, 808, 129, 366, 208, 385)),]

# filter for low and high autocor\
figDat_temp[figDat_temp$missingness=="MAR" & figDat_temp$autoCor <=0.3, "missingness"] <- "MCAR: Low AC"
figDat_temp[figDat_temp$missingness=="MAR" & figDat_temp$autoCor > 0.3 & figDat_temp$autoCor < 0.6, "missingness"] <- "MCAR: Med. AC"
figDat_temp[figDat_temp$missingness=="MAR" & figDat_temp$autoCor >=  0.6, "missingness"] <- "MCAR: High AC"

# reorder factor levels for plotting ##
figDat_temp <- figDat_temp %>% 
  mutate(type=fct_relevel(type,c("Data Deletion Simple", "Data Deletion CC" ,"Kalman filter", "Multiple imputations","brms")))

# figDat_temp %>% 
#   filter(missingprop_autocor != "y_noMiss",
#          param != "sigma") %>% 
#   filter(param == "phi") %>% 
#   ggplot() + 
#   facet_wrap(.~missingness) +
#   geom_point(aes(x = amtMiss, y = paramDiff, col = type), alpha = .2) + 
#   geom_smooth(aes(x = amtMiss, y = paramDiff, col = type)) + 
#   ylim(c(-5, 5))
#
## read in data for brms model runs w/ the complete, no-missing data simulations
figDat_brmsNoMiss <- read.csv("./data/model_results/gauss_sim_minMax_modelResults/AllParams_brmsNoMiss.csv") %>% 
  rename(simName = run_no,
         param = parameter, 
         param_value = mean, 
         param_se = sd,
         '2.5%' = X2.5.,
         '50%' = X50. ,
         '97.5%' = X97.5.) %>% 
  mutate(curSim = simName,
         autoCor = 0, 
         amtMiss = 0,
         param = replace(param, param %in% c("b_Intercept", "Intercept"), "intercept"),
         param = replace(param, param %in% c("b_light"), "light"),
         param = replace(param, param %in% c("b_discharge"), "discharge"),
  ) %>% 
  select(simName, curSim, missingprop_autocor, missingness, type, param, param_value, param_se, `2.5%`, `50%`, `97.5%`, autoCor, amtMiss) %>% 
  left_join(figDat_temp %>% select(simName, param, param_simVal) %>% unique()) %>% 
  mutate(paramDiff = ((param_value - param_simVal)/abs(param_simVal)),
         paramDiff_absDiff = abs((param_value - param_simVal)/abs(param_simVal)))

# because the model parameters are the same for all arima methods for time series w/ no missingness, then use the values from simple data deletion
# get figDat w/out y_noMiss data
figDat_missOnly <- figDat_temp %>% 
  filter(missingprop_autocor != "y_noMiss")
# get figDat with only y_noMiss data and remove duplicate model fits
figDat_noMiss <- figDat_temp %>% 
  filter(missingprop_autocor == "y_noMiss") %>% 
  filter(curSim < 1001)
figDat_noMiss <- figDat_noMiss %>% 
  rbind(figDat_noMiss %>% 
          filter(missingprop_autocor == "y_noMiss") %>% 
          mutate(type = replace(type, type == "dropNA_simple", "dropNA_complete")),
        figDat_noMiss %>% 
          filter(missingprop_autocor == "y_noMiss") %>% 
          mutate(type = replace(type, type == "dropNA_simple", "Kalman Filter")),
        figDat_noMiss %>% 
          filter(missingprop_autocor == "y_noMiss") %>% 
          mutate(type = replace(type, type == "dropNA_simple", "Multiple Imputations")),
        figDat_brmsNoMiss %>% 
          mutate(missingness = "MCAR: Med. AC")
  )
figDat_noMiss <- figDat_noMiss %>% 
  rbind(figDat_noMiss %>% 
          mutate(missingness = "MNAR"))

## add back together
figDat_temp <- figDat_missOnly %>% 
  rbind(figDat_noMiss)
# for the y_noMiss version, change amount missing to 0
figDat <- figDat_temp %>% 
  mutate(amtMiss = replace(amtMiss, missingprop_autocor == "y_noMiss", 0),
         missingness = replace(missingness,  missingprop_autocor == "y_noMiss" &  missingness =="MCAR: Low AC", "MCAR: Med. AC")) %>% 
  filter(param != "sigma",
         # missingprop_autocor != "y_noMiss"
  ) %>%
  # try removing simulations that have a parameter that is very, very small
  # filter( param_simVal > 0.05) %>% 
  mutate(autoCor = round(autoCor, 1), 
         #amtMiss = round(amtMiss,1),
         amtMiss = replace(amtMiss, amtMiss <=0.3 & amtMiss > 0, 0.2),
         amtMiss = replace(amtMiss, amtMiss > 0.3 & amtMiss <=0.5, 0.4),
         amtMiss = replace(amtMiss, amtMiss > 0.5, 0.6),
         param = replace(param, param %in% c("discharge", "light"), "Beta covariates"),
         param = replace(param, param == "intercept", "Intercept"),
         param = replace(param, param == "phi", "Phi")
  ) %>% 
  group_by(missingness, type, param, autoCor, amtMiss) %>% 
  summarize(paramDiffAbsDiff_mean = mean(paramDiff_absDiff, na.rm = TRUE),
            paramDiffAbsDiff_med = median(paramDiff_absDiff, na.rm = TRUE),
            paramDiffAbsDiff_SD = sd(paramDiff_absDiff, na.rm = TRUE),
            n_paramDiffAbsDiff = length(paramDiff_absDiff),
            paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n_paramDiff = length(paramDiff),
            SE_mean = mean(param_se, na.rm = TRUE) # the mean of the parameter standard error (not standardized, but maybe should be?)
  ) %>% 
  mutate(tooBigSD = ifelse(paramDiff_SD > 1.96, yes = 1, no = NA)) %>% 
  #filter(n_paramDiff  > 25)  %>% # drop combinations that have fewer than 300 observations
  filter(amtMiss <=.65, 
         param != "Intercept")

##make heat map ##
# Make heatmaps for Gaussian MCAR data - simulations -------------------------------------
# bin amt missing and autocorr (average paramDiff)
# figDat <- gauss_sim_figDat %>% 
#   filter(param != "sigma" & 
#            param != "intercept" &
#            missingness == "MAR") %>% 
#   mutate(
#          autoCor = round(autoCor, 1), 
#          #amtMiss = round(amtMiss, 1),
#         # no longer have continuously-spaced amt. missingness, so need to bin that rather than just rounding (can keep the autoCor the same, since those values were generated in an evenly-spaced continuum)
#         amtMiss = replace(amtMiss, amtMiss <=0.3, 0.2), 
#         amtMiss = replace(amtMiss, amtMiss > 0.3 & amtMiss <=0.5, 0.4), 
#         amtMiss = replace(amtMiss, amtMiss > 0.5, 0.6), 
#         param = replace(param, param %in% c("discharge", "light"), "Beta covariates"), # group light and discharge
#          param = replace(param, param == "intercept", "Intercept"),
#          param = replace(param, param == "phi", "Phi")) %>% 
#   group_by(missingness, type, param, autoCor, amtMiss) %>% 
#   summarize(paramDiffAbsDiff_mean = mean(paramDiff_absDiff, na.rm = TRUE),
#             paramDiffAbsDiff_med = median(paramDiff_absDiff, na.rm = TRUE),
#             paramDiffAbsDiff_SD = sd(paramDiff_absDiff, na.rm = TRUE),
#             n_paramDiffAbsDiff = length(paramDiff_absDiff), 
#             paramDiff_mean = mean(paramDiff, na.rm = TRUE),
#             paramDiff_med = median(paramDiff, na.rm = TRUE),
#             paramDiff_SD = sd(paramDiff, na.rm = TRUE),
#             n = length(paramDiff)) %>% 
#   #filter(n  > 100) %>% # drop combinations that have fewer than 100 observations
#   mutate(tooBigSD = ifelse(paramDiff_SD > 1.96, yes = 1, no = NA)) %>% 
#   filter(amtMiss <=.65)
# # only consider missingness of 50% or less

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
         param = replace(param, param == "Beta covariates", 'beta')) %>% 
  filter(missingness != "MNAR")

## make heatmap for median of parameter recovery
(heatMap_median_MCAR <- ggplot(data = figDat2, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiff_med), size=5) + 
    #scale_fill_viridis_c(name = "value" , option = "A") +
    scale_fill_gradient2(high = "darkblue",
                         mid = "orange",
                         low = "yellow" ,
                         midpoint = -.5, limits = c( min(figDat2$paramDiff_med),  max(figDat2$paramDiff_med)),#2.1),  
                         na.value = "lightgrey") +
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
    ylim(-0.1,1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Median error of parameter recovery")) 

## make heatmap for SE of parameter recovery
(heatMap_SE_MCAR <-ggplot(data = figDat2, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiffAbsDiff_med), size=5) + 
    #scale_fill_viridis_c(name = "value" , option = "D") +
    scale_fill_gradient2(low= "darkblue",
                         mid = "turquoise",
                         high = "green" ,
                         midpoint = .5, limits = c(0, max(figDat2$paramDiffAbsDiff_med)), 
                         na.value = "lightgrey") +
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
    ylim(-0.1,1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Median absolute error of parameter recovery")) 

## make heatmap for SD of parameter recovery
# calculate coverage 
figDat_cov_temp <- figDat_temp  %>% 
  mutate(CI95_lower = param_value - 1.96*param_se,
         CI95_upper = param_value + 1.96*param_se) %>% 
  filter(param != "sigma") %>% 
  filter(!is.na(param_se))# randomly there are some model runs that don't have SE? 

# is the true parameter within the 95% CI? 
figDat_cov_temp$coverage <- c(figDat_cov_temp$param_simVal >= figDat_cov_temp$CI95_lower & 
                                    figDat_cov_temp$param_simVal <= figDat_cov_temp$CI95_upper)

## count the # of models w/ and without coverage for each bin of missingness and autocorrelation
figDat_cov <- figDat_cov_temp %>% 
  mutate(amtMiss = replace(amtMiss, missingprop_autocor == "y_noMiss", 0),
         missingness = replace(missingness,  missingprop_autocor == "y_noMiss" &  missingness =="MCAR: Low AC", "MCAR: Med. AC")) %>% 
  filter(param != "sigma",
         # missingprop_autocor != "y_noMiss"
  ) %>%
  filter(param != "sigma",
         param != "intercept",
         amtMiss <=.65,) %>% 
  mutate(autoCor =  round(autoCor, 1), 
         #amtMiss = round(amtMiss,1),
         amtMiss = replace(amtMiss, amtMiss <=0.3 & amtMiss > 0, 0.2),
         amtMiss = replace(amtMiss, amtMiss > 0.3 & amtMiss <=0.5, 0.4),
         amtMiss = replace(amtMiss, amtMiss > 0.5, 0.6),
         param = replace(param, param %in% c("discharge", "light"), "beta"),
         param = replace(param, param == "intercept", "Intercept"),
         param = replace(param, param == "phi", "phi")
  ) %>% 
  group_by(missingness, type, param, autoCor, amtMiss) %>% 
  summarize(coverageNumber = sum(coverage), # the number of models that have coverage
            modelRunN = length(!is.na(coverage))# the total number of models 
            ) %>% 
  mutate(coveragePerc = coverageNumber/modelRunN)
    
figDat_covMAR <- figDat_cov %>% 
  filter(missingness != "MNAR")
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
   ylim(-0.1,1) +
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



# Get figure of MCAR median error for phi  simulations ----------------------------------
figDat3 <- figDat2 %>% 
  filter(param == "phi")
(heatMap_median_MCAR_slice <-ggplot(data = figDat3, aes(x=amtMiss, y=autoCor)) + 
   geom_tile(aes(fill=paramDiff_med), size=5) + 
   #scale_fill_viridis_c(name = "value" , option = "A") +
    scale_fill_gradient2(high = "darkblue",
                         mid = "orange",
                         low = "yellow" ,
                         midpoint = -.5, limits = c(min(figDat2$paramDiff_med), max(figDat2$paramDiff_med)),  na.value = "lightgrey") + # use same color scale as previously so figures are comprable
   #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
   ggh4x::facet_grid2(
     ~factor(type, levels = c("dropNA_complete", "dropNA_simple", "Kalman Filter", "Multiple Imputations", "brms"), 
             labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
     , labeller = label_parsed) +
   xlab("Proportion of missing data")+
   ylab("Autocorrelation in missingness") +
   guides(fill = guide_colorbar("Median Error")) +
   theme_classic() +
   # ylim(0,1) +
   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
   ggtitle("Missing completely at random: Median error of parameter recovery, \u03D5" )) 


png(file = "./figures/heatmap_GaussianMCAR_justPhi.png", width = 9, height = 3, units = "in", res = 700)
heatMap_median_MCAR_slice
dev.off()

# # Make heatmaps for Gaussian MNAR data  simulations ------------------------------------
# # bin amt missing and autocorr (average paramDiff)
# figDat <- gauss_sim_figDat %>% 
#   filter(param != "sigma" & 
#            param != "intercept" & 
#            missingness == "MNAR") %>% 
#   mutate(autoCor = 0, 
#          #amtMiss = round(amtMiss, 1),
#          # no longer have continuously-spaced amt. missingness, so need to bin that rather than just rounding (can keep the autoCor the same, since those values were generated in an evenly-spaced continuum)
#          amtMiss = replace(amtMiss, amtMiss <=0.3, 0.2), 
#          amtMiss = replace(amtMiss, amtMiss > 0.3 & amtMiss <=0.5, 0.4), 
#          amtMiss = replace(amtMiss, amtMiss > 0.5, 0.6), 
#          param = replace(param, param %in% c("discharge", "light"), "beta"), # group light and discharge
#          param = replace(param, param == "intercept", "Intercept"),
#          param = replace(param, param == "phi", "phi")) %>% 
#   group_by(missingness, type, param, autoCor, amtMiss) %>% 
#   summarize(paramDiffAbsDiff_mean = mean(paramDiff_absDiff, na.rm = TRUE),
#             paramDiffAbsDiff_med = median(paramDiff_absDiff, na.rm = TRUE),
#             paramDiffAbsDiff_SD = sd(paramDiff_absDiff, na.rm = TRUE),
#             n_paramDiffAbsDiff = length(paramDiff_absDiff),
#     paramDiff_mean = mean(paramDiff, na.rm = TRUE),
#             paramDiff_med = median(paramDiff, na.rm = TRUE),
#             paramDiff_SD = sd(paramDiff, na.rm = TRUE),
#             n = length(paramDiff)) %>% 
#   filter(n  > 300) %>% # drop combinations that have fewer than 300 observations
#   mutate(tooBigSD = ifelse(paramDiff_SD > 1.96, yes = 1, no = NA)) 
# 
# ## make heatmap for mean of parameter recovery
# (heatMap_medians_MNAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
#     geom_tile(aes(fill=paramDiff_med), size=5) + 
#     #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
#     ggh4x::facet_grid2(factor(param, levels = c("phi", "beta")
#     )
#     ~factor(type, levels = c("dropNA_complete", "dropNA_simple", "Kalman Filter", "Multiple Imputations", "brms"), 
#             labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
#     , labeller = label_parsed) +
#     scale_fill_viridis_c(begin=1, end=0, option = "A")+
#     xlab("Proportion of missing data")+
#     ylab("Autocorrelation in missingness") +
#     guides(fill = guide_colorbar("Median Error")) +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           axis.title.y = element_blank()) +
#     ggtitle("Median of parameter error across sims.")) 
# 
# ## make heatmap for mean of parameter recovery
# (heatMap_SE_MNAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
#     geom_tile(aes(fill=paramDiffAbsDiff_med), size=5) + 
#     #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
#     ggh4x::facet_grid2(factor(param, levels = c("phi", "beta")
#     )
#     ~factor(type, levels = c("dropNA_complete", "dropNA_simple", "Kalman Filter", "Multiple Imputations", "brms"), 
#             labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
#     , labeller = label_parsed) +
#     scale_fill_viridis_c(begin=1, end=0, option = "D")+
#     xlab("Proportion of missing data")+
#     ylab("Autocorrelation in missingness") +
#     theme_classic() +
#     guides(fill = guide_colorbar("Median \n Absolute Error")) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           axis.title.y = element_blank()) +
#     ggtitle("Absolute standard error of parameter recovery")) 
# 
# 
# figDat_covMNAR$autoCor <- 0
# (heatMap_cov_MNAR <- ggplot(data = figDat_covMNAR, aes(x=amtMiss, y=autoCor)) + 
#     geom_tile(aes(fill=coveragePerc*100), size=5) + 
#     ggh4x::facet_grid2(factor(param, levels = c("phi", "beta")
#     )
#     ~factor(type, levels = c("dropNA_complete", "dropNA_simple", "Kalman Filter", "Multiple Imputations", "brms"), 
#             labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
#     , labeller = label_parsed) +
#     #scale_fill_distiller(palette = "Greys", direction = 1, name = "value") +#value = SD
#     #scale_fill_viridis_c(begin=0, end=1, option = "G", name = "value" )+ 
#     scale_fill_gradient(low = "darkblue",
#                         high = "lightblue" ,
#                         limits = c(0,100),  na.value = "lightgrey") +
#     xlab("Proportion of missing data")+
#     guides(fill = guide_colorbar("% Coverage")) +
#     ylab("Autocorrelation in missingness") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     ggtitle("% of model runs where the 95% CI \n includes the simulation parameter"))
# 
# 
# ## save figures
# png(file = "./figures/heatmap_GaussianMNAR_median.png", width = 8, height = 6, units = "in", res = 700)
# heatMap_medians_MNAR
# dev.off()
# 
# png(file = "./figures/heatmap_GaussianMNAR_SE.png", width = 8, height = 6, units = "in", res = 700)
# heatMap_SE_MNAR
# dev.off()
# 
# png(file = "./figures/heatmap_GaussianMNAR_coverage.png", width = 8, height = 6, units = "in", res = 700)
# heatMap_cov_MNAR
# dev.off()
# 
# png(file = "./figures/heatmap_GaussianMNAR_all.png", width = 9, height = 12, units = "in", res = 700)
# ggarrange(heatMap_medians_MNAR, heatMap_SE_MNAR, heatMap_cov_MNAR, ncol = 1)
# dev.off()
# 
# 


# Make heatmaps for Gaussian MCAR data - empirical -------------------------------------
# # get MCAR data 
# gauss_real_arima <- read.csv("./data/model_results/gauss_real_MAR_arima_modResults/au_sable/gauss_auSable_real_MAR_arima_FORECASTvals.csv") %>% 
#   select(missingprop_autocor:run_no)
# gauss_real_brms <- read.csv("./data/model_results/gauss_real_MAR_brms_modResults/auSable/gauss_auSable_real_MAR_brms_FORECASTvals.csv") %>% 
#   select(missingprop_autocor:sd, missingness:run_no) %>% 
#   rename(parameters = parameter, 
#          param_value = mean, 
#          param_se = sd)
# 
# gauss_real_figDat <- gauss_real_arima %>% 
#   rbind(gauss_real_brms) %>% 
#   mutate(amtMiss = as.numeric(str_extract(missingprop_autocor, pattern = "[0-9/.0-9]+")),
#          autoCor = as.numeric(str_extract(missingprop_autocor, pattern = "[0-9/.0-9]+$")),
#          parameters = replace(parameters, parameters == "xreg1", "light"),
#          parameters = replace(parameters, parameters == "b_light","light"),
#          parameters = replace(parameters, parameters == "b_Intercept","intercept"),
#          parameters = replace(parameters, parameters == "b_Q","Q"),
#          parameters = replace(parameters, parameters == "xreg2","Q"),
#          parameters = replace(parameters, parameters == "ar1","phi")
#          )
# 
# # get MNAR data 
# gauss_real_arima_MNAR <- read.csv("./data/model_results/gauss_real_MNAR_arima_modResults/au_sable/gauss_auSable_real_MNAR_arima_FORECASTvals.csv") %>% 
#   select(missingprop_autocor:run_no)
# gauss_real_brms_MNAR <- read.csv("./data/model_results/gauss_real_MNAR_brms_modResults/auSable/brmsvals.csv") %>% 
#   select(missingprop_autocor:sd, missingness:run_no) %>% 
#   rename(parameters = parameter, 
#          param_value = mean, 
#          param_se = sd)
# 
# gauss_real_figDat_MNAR <- gauss_real_arima_MNAR %>% 
#   rbind(gauss_real_brms_MNAR) %>% 
#   mutate(autoCor = NA,
#     amtMiss = as.numeric(str_extract(missingprop_autocor, pattern = "[0-9/.0-9]+$"))
#   ) %>% 
#   mutate(amtMiss = replace(amtMiss, is.na(amtMiss), 0.44)) %>% 
#   mutate(
#          parameters = replace(parameters, parameters == "xreg1", "light"),
#          parameters = replace(parameters, parameters == "b_light","light"),
#          parameters = replace(parameters, parameters == "b_Intercept","intercept"),
#          parameters = replace(parameters, parameters == "Intercept","intercept"),
#          parameters = replace(parameters, parameters == "b_discharge","Q"),
#          parameters = replace(parameters, parameters == "xreg2","Q"),
#          parameters = replace(parameters, parameters == "ar1","phi")
#   )
# ## add the MNAR data to the MAR data
# gauss_real_figDat <- gauss_real_figDat %>% 
#   rbind(gauss_real_figDat_MNAR)
# ## get the data from models fit to the 'complete' dataset (has a few NAs), which will be the parameters that those fit to the missing datasets will be compared to
# gauss_real_RefParams_arima <- read.csv("./data/model_results/gauss_real_MAR_arima_modResults/au_sable/CompleteDataset_arimavals.csv") %>% 
#   select(missingprop_autocor:run_no)
# gauss_real_RefParams_brms <- read.csv("./data/model_results/gauss_real_MAR_brms_modResults/auSable/CompleteDataset_brmsvals.csv") %>% 
#   select(missingprop_autocor:sd, missingness:run_no) %>% 
#   rename(parameters = parameter, 
#          param_value = mean, 
#          param_se = sd) %>% 
#   filter(parameters != "Intercept")
# 
# gauss_real_RefParams <- gauss_real_RefParams_arima %>% 
#   rbind(gauss_real_RefParams_brms) %>% 
#   mutate(
#     parameters = replace(parameters, parameters == "xreg1", "light"),
#     parameters = replace(parameters, parameters == "b_light","light"),
#     parameters = replace(parameters, parameters == "b_Intercept","intercept"),
#     parameters = replace(parameters, parameters == "b_Q","Q"),
#     parameters = replace(parameters, parameters == "xreg2","Q"),
#     parameters = replace(parameters, parameters == "ar1","phi")
#   ) %>% 
#   rename(param_simVal = param_value)
# 
# ## add the parameters from the models fit to the 'complete' dataset as reference 
# gauss_real_figDat <- gauss_real_figDat %>% 
#   left_join(gauss_real_RefParams %>% select(parameters, param_simVal, type)) %>% 
#   # calculate difference between model-derived parameters w/ missingness and no missingness
#   dplyr::mutate(paramDiff = ((param_value - param_simVal)/abs(param_simVal)),
#                 paramDiff_absDiff = abs((param_value - param_simVal)/abs(param_simVal))) 
# ## save for later 
# saveRDS(gauss_real_figDat, "./data/model_results/gauss_real_ModelResults.rds")

# # ggplot(gauss_real_figDat) + 
# #   geom_boxplot(aes(x = parameters, y = param_value, col = type))
# # bin amt missing and autocorr (average paramDiff)
# figDat <- gauss_real_figDat %>% 
#   filter(parameters != "sigma" & 
#            parameters != "intercept" &
#            parameters != "Intercept" &
#            missingness == "MAR") %>% 
#   mutate(
#     autoCor = round(autoCor, 1), 
#     #amtMiss = round(amtMiss, 1),
#     # no longer have continuously-spaced amt. missingness, so need to bin that rather than just rounding (can keep the autoCor the same, since those values were generated in an evenly-spaced continuum)
#     amtMiss = replace(amtMiss, amtMiss <=0.3, 0.2), 
#     amtMiss = replace(amtMiss, amtMiss > 0.3 & amtMiss <=0.5, 0.4), 
#     amtMiss = replace(amtMiss, amtMiss > 0.5, 0.6), 
#     parameters = replace(parameters, parameters %in% c("discharge", "light", "Q"), "Beta covariates"), # group light and discharge
#     parameters = replace(parameters, parameters == "intercept", "Intercept"),
#     parameters = replace(parameters, parameters == "phi", "Phi")) %>% 
#   dplyr::mutate(paramDiff = ((param_value - param_simVal)/abs(param_simVal)),
#                 paramDiff_absDiff = abs((param_value - param_simVal)/abs(param_simVal))) %>% 
#   group_by(missingness, type, parameters, autoCor, amtMiss) %>% 
#   summarize(paramDiffAbsDiff_mean = mean(paramDiff_absDiff, na.rm = TRUE),
#             paramDiffAbsDiff_med = median(paramDiff_absDiff, na.rm = TRUE),
#             paramDiffAbsDiff_SD = sd(paramDiff_absDiff, na.rm = TRUE),
#             n_paramDiffAbsDiff = length(paramDiff_absDiff), 
#             paramDiff_mean = mean(paramDiff, na.rm = TRUE),
#             paramDiff_med = median(paramDiff, na.rm = TRUE),
#             paramDiff_SD = sd(paramDiff, na.rm = TRUE),
#             n = length(paramDiff)) %>% 
#   #filter(n  > 100) %>% # drop combinations that have fewer than 100 observations
#   mutate(tooBigSD = ifelse(paramDiff_SD > 1.96, yes = 1, no = NA)) %>% 
#   filter(amtMiss <=.65)
# # only consider missingness of 50% or less
# 
# # update names for missing data approach
# type_names <- c(
#   "dropNA_complete" = "Data Deletion-Complete",
#   "dropNA_simple" = "Data Deletion-Simple",
#   "Kalman Filter" =  "Kalman filter",
#   "Multiple Imputations" = "Multiple Imputation",
#   "brms" = "Data Augmentation", 
#   "Phi" = "Phi", 
#   "Beta covariates" = "Beta covariates"
# )
# 
# figDat2 <- figDat %>% 
#   mutate(parameters = replace(parameters, parameters == "Phi", "phi"), 
#          parameters = replace(parameters, parameters == "Beta covariates", 'beta'))
# 
# ## make heatmap for median of parameter recovery
# (heatMap_median_MCAR <-ggplot(data = figDat2, aes(x=amtMiss, y=autoCor)) + 
#     geom_tile(aes(fill=paramDiff_med), size=5) + 
#     #scale_fill_viridis_c(name = "value" , option = "A") +
#     scale_fill_gradient2(high = "darkblue",
#                          mid = "orange",
#                          low = "yellow" ,
#                          midpoint = -.5, limits = c( min(figDat2$paramDiff_med),  max(figDat2$paramDiff_med)),#2.1),  
#                          na.value = "lightgrey") +
#     #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
#     ggh4x::facet_grid2(factor(parameters, levels = c("phi", "beta")
#     )
#     ~factor(type, levels = c("dropNA_complete", "dropNA_simple", "Kalman Filter", "Multiple Imputations", "brms"), 
#             labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
#     , labeller = label_parsed) +
#     xlab("Proportion of missing data")+
#     ylab("Autocorrelation in missingness") +
#     guides(fill = guide_colorbar("Median Error")) +
#     theme_classic() +
#     ylim(0,1) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     ggtitle("Median error of parameter recovery")) 
# 
# ## make heatmap for SE of parameter recovery
# (heatMap_SE_MCAR <-ggplot(data = figDat2, aes(x=amtMiss, y=autoCor)) + 
#     geom_tile(aes(fill=paramDiffAbsDiff_med), size=5) + 
#     #scale_fill_viridis_c(name = "value" , option = "D") +
#     scale_fill_gradient2(low= "darkblue",
#                          mid = "turquoise",
#                          high = "green" ,
#                          midpoint = .5, limits = c(0, max(figDat2$paramDiffAbsDiff_med)), 
#                          na.value = "lightgrey") +
#     #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
#     ggh4x::facet_grid2(factor(parameters, levels = c("phi", "beta")
#     )
#     ~factor(type, levels = c("dropNA_complete", "dropNA_simple", "Kalman Filter", "Multiple Imputations", "brms"), 
#             labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
#     , labeller = label_parsed) +
#     xlab("Proportion of missing data")+
#     ylab("Autocorrelation in missingness") +
#     guides(fill = guide_colorbar("Median \nAbsolute Error")) +
#     theme_classic() +
#     ylim(0,1) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     ggtitle("Median absolute error of parameter recovery")) 
# 
# ## make heatmap for SD of parameter recovery
# # calculate coverage 
# figDat_cov_temp <- gauss_real_figDat  %>% 
#   mutate(CI95_lower = param_value - 1.96*param_se,
#          CI95_upper = param_value + 1.96*param_se) %>% 
#   filter(parameters != "sigma") %>% 
#   filter(!is.na(param_se))# randomly there are some model runs that don't have SE? 
# 
# # is the true parameter within the 95% CI? 
# figDat_cov_temp$coverage <- c(figDat_cov_temp$param_simVal >= figDat_cov_temp$CI95_lower & 
#                                 figDat_cov_temp$param_simVal <= figDat_cov_temp$CI95_upper)
# 
# ## count the # of models w/ and without coverage for each bin of missingness and autocorrelation
# figDat_cov <- figDat_cov_temp %>% 
#   filter(parameters != "sigma",
#          parameters != "intercept",
#          amtMiss <=.5,) %>% 
#   mutate(autoCor = round(autoCor, 1), 
#          #amtMiss = round(amtMiss, 1), # no longer have continuously-spaced amt. missingness, so need to bin that rather than just rounding (can keep the autoCor the same, since those values were generated in an evenly-spaced continuum)
#          amtMiss = replace(amtMiss, amtMiss <=0.3, 0.2), 
#          amtMiss = replace(amtMiss, amtMiss > 0.3 & amtMiss <=0.5, 0.4), 
#          amtMiss = replace(amtMiss, amtMiss > 0.5, 0.6), 
#          parameters = replace(parameters, parameters %in% c("discharge", "light"), "beta"),
#          #parameters = replace(parameters, parameters == "intercept", "Intercept"),
#          parameters = replace(parameters, parameters == "phi", "phi")
#   ) %>% 
#   group_by(missingness, type, parameters, autoCor, amtMiss) %>% 
#   summarize(coverageNumber = sum(coverage), # the number of models that have coverage
#             modelRunN = length(!is.na(coverage))# the total number of models 
#   ) %>% 
#   mutate(coveragePerc = coverageNumber/modelRunN)
# 
# figDat_covMAR <- figDat_cov %>% 
#   filter(missingness == "MAR")
# figDat_covMNAR <- figDat_cov %>% 
#   filter(missingness %in% c("MNAR"))
# 
# (heatMap_cov_MCAR <- ggplot(data = figDat_covMAR, aes(x=amtMiss, y=autoCor)) + 
#     geom_tile(aes(fill=coveragePerc*100), size=5) + 
#     #geom_tile(data = figDat_covMAR[figDat_covMAR$coveragePerc >.96,], aes(), fill = "grey") + 
#     ggh4x::facet_grid2(factor(parameters, levels = c("phi", "beta")
#     )
#     ~factor(type, levels = c("dropNA_complete", "dropNA_simple", "Kalman Filter", "Multiple Imputations", "brms"), 
#             labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
#     , labeller = label_parsed) +
#     #scale_fill_distiller(palette = "Greys", direction = 1, name = "value") +#value = SD
#     #scale_fill_viridis_c(begin=0, end=1, option = "G", name = "value" )+ 
#     scale_fill_gradient(low = "darkblue",
#                         high = "lightblue" ,
#                         limits = c(0,100),  na.value = "lightgrey") +
#     guides(fill = guide_colorbar("Median \n Absolute Error")) +
#     xlab("Proportion of missing data")+
#     ylab("Autocorrelation in missingness") +
#     guides(fill = guide_colorbar("% Coverage")) +
#     theme_classic() +
#     ylim(0,1) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     ggtitle("% of model runs where the 95% CI includes the simulation parameter"))
# 
# ## save figures
# png(file = "./figures/heatmap_GaussianMCAR_median.png", width = 8, height = 6, units = "in", res = 700)
# heatMap_median_MCAR
# dev.off()
# 
# png(file = "./figures/heatmap_GaussianMCAR_SE.png", width = 8, height = 6, units = "in", res = 700)
# heatMap_SE_MCAR
# dev.off()
# 
# png(file = "./figures/heatmap_GaussianMCAR_coverage.png", width = 8, height = 6, units = "in", res = 700)
# heatMap_cov_MCAR
# dev.off()
# 
# png(file = "./figures/heatmap_GaussianMCAR_all.png", width = 9.5, height = 12, units = "in", res = 700)
# ggarrange(heatMap_median_MCAR, heatMap_SE_MCAR, heatMap_cov_MCAR, ncol = 1) %>% 
#   annotate_figure(top = text_grob("Parameter recovery from real-valued time series with MCAR data", just = .65, 
#                                   size = 16, face = "bold"))
# dev.off()
# 
# 
# 
# # Get figure of MCAR median error for phi -  empirical ----------------------------------
# figDat3 <- figDat2 %>% 
#   filter(parameters == "phi")
# (heatMap_median_MCAR_slice <-ggplot(data = figDat3, aes(x=amtMiss, y=autoCor)) + 
#     geom_tile(aes(fill=paramDiff_med), size=5) + 
#     #scale_fill_viridis_c(name = "value" , option = "A") +
#     scale_fill_gradient2(high = "darkblue",
#                          mid = "orange",
#                          low = "yellow" ,
#                          midpoint = -.5, limits = c(min(figDat2$paramDiff_med), max(figDat2$paramDiff_med)),  na.value = "lightgrey") + # use same color scale as previously so figures are comprable
#     #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
#     ggh4x::facet_grid2(
#       ~factor(type, levels = c("dropNA_complete", "dropNA_simple", "Kalman Filter", "Multiple Imputations", "brms"), 
#               labels = c('"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Kalman filter"', '"Multiple Imputation"', '"Data Augmentation"'))
#       , labeller = label_parsed) +
#     xlab("Proportion of missing data")+
#     ylab("Autocorrelation in missingness") +
#     guides(fill = guide_colorbar("Median Error")) +
#     theme_classic() +
#     ylim(0,1) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     ggtitle("Missing completely at random: Median error of parameter recovery, \u03D5" )) 
# 
# 
# png(file = "./figures/heatmap_GaussianMCAR_justPhi.png", width = 9, height = 3, units = "in", res = 700)
# heatMap_median_MCAR_slice
# dev.off()
