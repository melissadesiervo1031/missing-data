## draft of the Population heat map figure ##

# Load packages
library(tidyverse)
library(ggpubr)
library(ggh4x)

## read in data (prepared in "Population_output_plots.R" script)

ricDat_new <- readRDS( "./data/model_results/ricker_sim_ParamsForFigures.rds")

ricDat_new$rowID <- seq(1:nrow(ricDat_new))


### prepare data 

#make into long data.frame 
paramEstLong <- ricDat_new %>% 
  #select(SimNumber, r_sim, actAutoCorr, actPropMiss, type, r_est, r_se, status, r_paramDiff) %>% 
  pivot_longer(cols = c(r_est, alpha_est), 
               values_to = "paramEst", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  mutate(param = str_replace(param, "_est", "")) %>% 
  select(-r_sim, -alpha_sim, -N0_sim, -r_se, -alpha_se)

paramSimLong <- ricDat_new %>% 
  pivot_longer(cols = c(r_sim, alpha_sim), 
               values_to = "paramSim", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  mutate(param = str_replace(param, "_sim", "")) %>% 
  select(-r_est, -alpha_est, -N0_sim, -r_se, -alpha_se) %>% 
  unique()

paramSELong <- ricDat_new %>% 
  pivot_longer(cols = c(r_se, alpha_se), 
               values_to = "paramSE", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  mutate(param = str_replace(param, "_se", "")) %>% 
  select(-r_sim, -alpha_sim, -N0_sim, -r_est, -alpha_est) %>% 
  # remove the values for Expectation Maximization, since we don't have SE for that method
  filter(type != "ExpectationMaximization")


ricDat_new_long <- left_join(paramEstLong, paramSimLong) %>% 
  left_join(paramSELong)

#calculate standardized parameter estimates
ricDat_new_long <- ricDat_new_long %>%
  mutate("paramDiff" = (paramEst - paramSim)/abs(paramSim),
         "paramDiff_abs" = abs(paramEst - paramSim)/abs(paramSim))

# # filter for low and high autocor
# ricDat_new_long[ricDat_new_long$autoCorr <=0.3 & !is.na(ricDat_new_long$autoCorr), "missingnessType"] <- "MCAR"
# ricDat_new_long[ricDat_new_long$autoCorr >0.3 & ricDat_new_long$autoCorr <0.6 & !is.na(ricDat_new_long$autoCorr), "missingnessType"] <- "MCAR"
# ricDat_new_long[ricDat_new_long$autoCorr  >= 0.6 & !is.na(ricDat_new_long$autoCorr), "missingnessType"] <- "MCAR"


ricDat_new_long$type <- factor(ricDat_new_long$type, levels = c("DataAugmentation", "CompleteCaseDropNA", 
                                                                "dropNA", "ExpectationMaximization", "MultipleImputations"))

# for data w/ no missingness, replace NA in propMiss with 0
ricDat_new_long[ricDat_new_long$listName == "y_noMiss", "propMiss"] <- 0
ricDat_new_long[ricDat_new_long$listName == "y_noMiss", "missingnessType"] <- "MCAR"
ricDat_new_long_NOMISS <- ricDat_new_long[ricDat_new_long$listName == "y_noMiss", ]
ricDat_new_long_NOMISS$missingnessType <- "MNAR"
ricDat_new_long <- ricDat_new_long %>% 
  rbind(ricDat_new_long_NOMISS)
ricDat_new_long_NOMISS_ForMI <-  ricDat_new_long[ricDat_new_long$listName == "y_noMiss", ]
ricDat_new_long_NOMISS_ForMI <- ricDat_new_long_NOMISS_ForMI %>% 
  filter(type == "dropNA") %>% 
  mutate(type = "MultipleImputations")
ricDat_new_long <- ricDat_new_long %>% 
  rbind(ricDat_new_long_NOMISS_ForMI)
# if the model run is MCAR and includes no missing data, list autCorr as 0
ricDat_new_long[is.na(ricDat_new_long$autoCorr) & ricDat_new_long$missingnessType == "MCAR", "autoCorr"] <- 0
# reformat data 
figDat_lines <- ricDat_new_long %>% 
  mutate(autoCor = round(autoCorr, 1), 
         amtMiss = propMiss,
         amtMiss = replace(amtMiss, propMiss <=0.3 & propMiss > 0, 0.2),
         amtMiss = replace(amtMiss, propMiss > 0.3 & propMiss <=0.5, 0.4),
         amtMiss = replace(amtMiss, propMiss > 0.5, 0.6)#,
         #amtMiss = round(propMiss, 1)
  ) %>% 
  group_by(missingnessType, type, param, amtMiss, autoCor) %>% 
  
  summarize(paramDiffAbsDiff_mean = mean(paramDiff_abs, na.rm = TRUE),
            paramDiffAbsDiff_med = median(paramDiff_abs, na.rm = TRUE),
            paramDiffAbsDiff_SD = sd(paramDiff_abs, na.rm = TRUE),
            n_paramDiffAbsDiff = length(paramDiff_abs),
            paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n_paramDiff = length(paramDiff),
            SE_mean = mean(paramSE, na.rm = TRUE) # the mean of the parameter standard error (not standardized, but maybe should be?)
  ) %>% 
  
  #filter(n  > 50)  %>% # drop combinations that have fewer than 300 observations
  filter(amtMiss <=.65)
# make types in the 'correct' order


# reorder factor levels for plotting ##
figDat_lines <- figDat_lines %>% 
  mutate(type=fct_relevel(type,c("dropNA", "CompleteCaseDropNA" ,"MultipleImputations","ExpectationMaximization","DataAugmentation")))


#figDat_lines2<-figDat_lines%>% filter(missingnessType %in% c("MCAR: Med. AC", "MNAR"))

### 

# #make into long data.frame 
# paramEstLong <- ricDat_new %>% 
#   #select(SimNumber, r_sim, actAutoCorr, actPropMiss, type, r_est, r_se, status, r_paramDiff) %>% 
#   pivot_longer(cols = c(r_est, alpha_est), 
#                values_to = "paramEst", 
#                names_to = "param", 
#                names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
#   select(-r_sim, -alpha_sim, -N0_sim, -r_se, -alpha_se)
# 
# paramSimLong <- ricDat_new %>% 
#   pivot_longer(cols = c(r_sim, alpha_sim), 
#                values_to = "paramSim", 
#                names_to = "param", 
#                names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
#   select(-r_est, -alpha_est, -N0_sim, -r_se, -alpha_se) %>% 
#   unique()
# 
# paramSELong <- ricDat_new %>% 
#   pivot_longer(cols = c(r_se, alpha_se), 
#                values_to = "paramSE", 
#                names_to = "param", 
#                names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
#   select(-r_sim, -alpha_sim, -N0_sim, -r_est, -alpha_est
#          ) %>% 
#   # remove the values for Expectation Maximization, since we don't have SE for that method
#   filter(type != "ExpectationMaximization")
# 
# 
# ricDat_new_long <- left_join(paramEstLong, paramSimLong) %>% 
#   left_join(paramSELong)
# 
# #calculate standardized parameter estimates
# ricDat_new_long <- ricDat_new_long %>%
#   mutate("paramDiff" = (paramEst - paramSim)/abs(paramSim),
#          "paramDiff_abs" = abs(paramEst - paramSim)/abs(paramSim))
# 
# # filter for low and high autocor
# ricDat_new_long[ricDat_new_long$autoCorr <=0.3 & !is.na(ricDat_new_long$autoCorr), "missingness"] <- "MCAR: Low AC"
# ricDat_new_long[ricDat_new_long$autoCorr >0.3 & ricDat_new_long$autoCorr <0.6 & !is.na(ricDat_new_long$autoCorr), "missingness"] <- "MCAR: Med. AC"
# ricDat_new_long[ricDat_new_long$autoCorr  >= 0.6 & !is.na(ricDat_new_long$autoCorr), "missingness"] <- "MCAR: High AC"
# 
# 
# ricDat_new_long$type <- factor(ricDat_new_long$type, levels = c("DataAugmentation", "CompleteCaseDropNA", 
#                                                                 "dropNA", "ExpectationMaximization", "MultipleImputations"))

figDat <- figDat_lines %>% 
  filter(missingnessType == "MCAR")
##make heat map! ##
# Make heatmaps for Gaussian MAR data -------------------------------------
# bin amt missing and autocorr (average paramDiff)
#%>%
  # mutate(autoCorr = round(autoCor, 1), 
  #        propMiss = round(propMiss, 1)) %>% 
  # group_by(missingness, type, param, autoCor, propMiss) %>% 
  # summarize(paramDiffAbsDiff_mean = mean(paramDiff_abs, na.rm = TRUE),
  #           paramDiffAbsDiff_med = median(paramDiff_abs, na.rm = TRUE),
  #           paramDiffAbsDiff_SD = sd(paramDiff_abs, na.rm = TRUE),
  #           n_paramDiffAbsDiff = length(paramDiff_abs),
  #           paramDiff_mean = mean(paramDiff, na.rm = TRUE),
  #           paramDiff_med = median(paramDiff, na.rm = TRUE),
  #           paramDiff_SD = sd(paramDiff, na.rm = TRUE),
  #           n = length(paramDiff)) %>% 
  # #filter(n  > 100) %>% # drop combinations that have fewer than 100 observations
  # filter(propMiss <=.5)
# only consider missingness of 50% or less


## make heatmap for mean of parameter recovery
(heatMap_median_MAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    facet_wrap(~missingnessType) + 
    geom_tile(aes(fill=paramDiff_med), size=5) + 
    #scale_fill_viridis_c(name = "value" , option = "A", direction = -1) +
    scale_fill_gradient2(low = "darkblue",
                         mid = "orange",
                         high = "yellow" ,
                         midpoint = .7, #limits = c(-.2, 1.3),#c(-.050, 0.82),  
                         na.value = "lightgrey") +
    #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    #facet_grid(~factor(figDat$param, levels = c("alpha", "r")) ~ type) +
    ggh4x::facet_grid2(factor(param, levels = c("alpha", "r")
    )
    ~factor(type, levels = 
            c("DataAugmentation", "CompleteCaseDropNA", "dropNA", "ExpectationMaximization", "MultipleImputations"),
            labels = c('"Data Augmentation"', '"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Expectation Maximization"', '"Multiple Imputation"' ))
    , labeller = label_parsed) +
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    guides(fill = guide_colorbar("Median Error")) +
    theme_classic() +
    #xlim(.05, .55) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Median error of parameter recovery")) 

## make heatmap for mean of *absolute value* of parameter recovery
(heatMap_SE_MAR <-ggplot(data = figDat, aes(x=amtMiss, y=autoCor)) + 
    geom_tile(aes(fill=paramDiffAbsDiff_med), size=5) + 
    #scale_fill_viridis_c(name = "value" , option = "D") +
    scale_fill_gradient2(low = "darkblue",
                         mid = "turquoise",
                        high = "green" ,
                        midpoint = .6, #limits = c(0, 0.82),  
                        na.value = "lightgrey") +
    #scale_fill_distiller(palette = "Spectral", direction = 1, name = "value")+
    #facet_grid(~factor(figDat$param, levels = c("alpha", "r")) ~ type) +
    ggh4x::facet_grid2(factor(param, levels = c("alpha", "r")
    )
    ~factor(type, levels = 
              c("DataAugmentation", "CompleteCaseDropNA", "dropNA", "ExpectationMaximization", "MultipleImputations"),
            labels = c('"Data Augmentation"', '"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Expectation Maximization"', '"Multiple Imputation"' ))
    , labeller = label_parsed) +
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    guides(fill = guide_colorbar("Median \nAbsolute Error")) +
    theme_classic() +
    #xlim(.05, .55) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Median absolute error of parameter recovery")) 

## make heatmap for coverage for parameter recovery
figDat_cov_temp <- ricDat_new_long %>% 
  mutate(CI95_lower = paramEst - 1.96*paramSE,
         CI95_upper = paramEst + 1.96*paramSE) %>% 
  filter(param != "sigma") %>% 
  filter(!is.na(paramSE))# randomly there are some model runs that don't have SE? 

# is the true parameter within the 95% CI? 
figDat_cov_temp$coverage <- c(figDat_cov_temp$paramSim >= figDat_cov_temp$CI95_lower & 
                                figDat_cov_temp$paramSim <= figDat_cov_temp$CI95_upper)

## count the # of models w/ and without coverage for each bin of missingness and autocorrelation
figDat_cov <- figDat_cov_temp %>% 
  filter(param != "sigma",
         param != "intercept",
         #propMiss <=.5
         ) %>% 
  mutate(autoCorr = round(autoCorr, 1), 
         amtMiss = propMiss,
         amtMiss = replace(amtMiss, propMiss <=0.3 & propMiss > 0, 0.2),
         amtMiss = replace(amtMiss, propMiss > 0.3 & propMiss <=0.5, 0.4),
         amtMiss = replace(amtMiss, propMiss > 0.5, 0.6)#,
  ) %>% 
  group_by(missingnessType, type, param, autoCorr, amtMiss) %>% 
  summarize(coverageNumber = sum(coverage), # the number of models that have coverage
            modelRunN = length(!is.na(coverage))# the total number of models 
  ) %>% 
  mutate(coveragePerc = coverageNumber/modelRunN)

(heatMap_cov_MAR <- ggplot(data = figDat_cov, aes(x=amtMiss, y=autoCorr)) + 
    geom_tile(aes(fill=coveragePerc*100), size=5) + 
    #facet_grid(~factor(figDat_cov$param, levels = c("alpha", "r")) ~ type) +
    ggh4x::facet_grid2(factor(param, levels = c("alpha", "r")
    )
    ~factor(type, levels = 
              c("DataAugmentation", "CompleteCaseDropNA", "dropNA", "ExpectationMaximization", "MultipleImputations"),
            labels = c('"Data Augmentation"', '"Data Deletion-Complete"', '"Data Deletion-Simple"', '"Expectation Maximization"', '"Multiple Imputation"' ))
    , labeller = label_parsed) +
    #scale_fill_distiller(palette = "Greys", direction = 1, name = "value") +
  scale_fill_gradient(low = "darkblue",
                       high = "lightblue" ,
                       limits = c(0,100),  na.value = "lightgrey") +
   # scale_fill_viridis_c(begin=0, end=1, option = "G", name = "value" )+
    xlab("Proportion of missing data")+
    ylab("Autocorrellation in missingness") +
    #xlim(.05, .55) +
    theme_classic() +
    guides(fill = guide_colorbar("% Coverage")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("% of model runs where the 95% CI includes the simulation parameter"))

## save figures
png(file = "./figures/heatmap_PoissonMCAR_median.png", width = 7, height = 6, units = "in", res = 700)
heatMap_median_MAR
dev.off()

png(file = "./figures/heatmap_PoissonMCAR_SE.png", width = 7, height = 6, units = "in", res = 700)
heatMap_SE_MAR
dev.off()

png(file = "./figures/heatmap_PoissonMCAR_coverage.png", width = 7, height = 6, units = "in", res = 700)
heatMap_cov_MAR
dev.off()

png(file = "./figures/heatmap_PoissonMCAR_all.png", width = 9.5, height = 12, units = "in", res = 700)
ggarrange(heatMap_median_MAR, heatMap_SE_MAR, ggarrange(heatMap_cov_MAR, ggplot() + geom_point(), nrow = 1, widths = c(.99, .01)), ncol = 1) %>% 
  annotate_figure(top = text_grob("Parameter recovery from time series of counts with MCAR data", just = .65, 
                                  size = 16, face = "bold"))
dev.off()

