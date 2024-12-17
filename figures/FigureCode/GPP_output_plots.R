#/////////////////
# This script makes figures for GPP simulated data 
# 30 October 2023
#/////////////////

# Load packages
library(tidyverse)
library(ggpubr)
library(forecast)
library(RColorBrewer)

## read in data 
figDat_temp <- readRDS("./data/model_results/gauss_sim_ModelResults.rds")

# remove data for simluation 376... has one really small parameter, which is causing a lot of outliers

figDat_temp <- figDat_temp[!(figDat_temp$simName %in% c(376, 831, 816, 461, 808, 129, 366, 208, 385)),]

# filter for low and high autocor\
figDat_temp[figDat_temp$missingness=="MAR" & figDat_temp$autoCor <=0.3, "missingness"] <- "MAR: Low AC"
figDat_temp[figDat_temp$missingness=="MAR" & figDat_temp$autoCor > 0.3 & figDat_temp$autoCor < 0.6, "missingness"] <- "MAR: Med. AC"
figDat_temp[figDat_temp$missingness=="MAR" & figDat_temp$autoCor >=  0.6, "missingness"] <- "MAR: High AC"

#
figDat_lines <- figDat_temp %>% 
  filter(param != "sigma") %>%
  # try removing simulations that have a parameter that is very, very small
  # filter( param_simVal > 0.05) %>% 
  mutate(autoCor = round(autoCor, 1), 
         amtMiss = round(amtMiss, 1),
         param = replace(param, param %in% c("discharge", "light"), "Beta covariates"),
         param = replace(param, param == "intercept", "Intercept"),
         param = replace(param, param == "phi", "Phi")
  ) %>% 
  group_by(missingness, type, param, amtMiss) %>% 
  summarize(paramDiffAbsDiff_mean = mean(paramDiff_absDiff, na.rm = TRUE),
            paramDiffAbsDiff_med = median(paramDiff_absDiff, na.rm = TRUE),
            paramDiffAbsDiff_SD = sd(paramDiff_absDiff, na.rm = TRUE),
            n_paramDiffAbsDiff = length(paramDiff_absDiff),
            paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n_paramDiff = length(paramDiff),
            SE_mean = mean(SE, na.rm = TRUE) # the mean of the parameter standard error (not standardized, but maybe should be?)
  ) %>% 
  
  #filter(n  > 50)  %>% # drop combinations that have fewer than 300 observations
  filter(amtMiss <=.5, 
         param != "Intercept")

figDat_linesMAR <- figDat_lines %>% 
  filter(missingness %in% c("MAR: High AC", "MAR: Low AC",  "MAR: Med. AC"))
figDat_linesMNAR <- figDat_lines %>% 
  filter(missingness %in% c("MNAR"))

# parameter recovery bias -------------------------------------------------
(gauss_paramRecovery_bias_MAR <- ggplot(data = figDat_linesMAR, aes(x = amtMiss, y = paramDiff_med)) +
   facet_grid(~factor(param, levels = c( "Intercept","Phi", "Beta covariates")) 
              ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC", "MNAR")), 
                scales = "free_y") + 
   geom_hline(aes(yintercept = 0), colour = "grey") + 
   #geom_errorbar(aes(ymin=paramDiff_mean - paramDiff_SD, ymax=paramDiff_mean + paramDiff_SD, color = as.factor(type)), 
   #size=0.3, width=0, position = position_dodge(width=0.03))+
   #geom_ribbon(aes(ymin = paramDiff_mean - paramDiff_SD, ymax = paramDiff_mean + paramDiff_SD, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
   geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
   #geom_point(data = data.frame("x" = c(0,.5), "y" = c(-.25, .1)), aes(x = x, y = y), alpha = 0.0000001) + ## add invisible points to change scales
   geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
   theme_classic() +
   xlab("Proportion of missing data")+ 
   theme(legend.position="top")+
   theme(legend.title=element_blank())+
   ylab("Median of parameter bias across sims.")+ 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
   xlim(c(-0.03,0.55)) + 
   #ylim(c(0,.3)) +
   scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#D95F02", "#7570B3",  "#E6AB02"),
                        labels = c("Data Augmentation", "Data Deletion-Complete","Data Deletion-Simple", "Kalman Filter", "Multiple Imputations"))
)
(gauss_paramRecovery_bias_MNAR <- ggplot(data = figDat_linesMNAR, aes(x = amtMiss, y = paramDiff_med)) +
    facet_grid(~factor(param, levels = c( "Intercept","Phi", "Beta covariates")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC", "MNAR")), 
               scales = "free_y") + 
    geom_hline(aes(yintercept = 0), colour = "grey") + 
    #geom_errorbar(aes(ymin=paramDiff_mean - paramDiff_SD, ymax=paramDiff_mean + paramDiff_SD, color = as.factor(type)), 
    #size=0.3, width=0, position = position_dodge(width=0.03))+
    #geom_ribbon(aes(ymin = paramDiff_mean - paramDiff_SD, ymax = paramDiff_mean + paramDiff_SD, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    #geom_point(data = data.frame("x" = c(0,.5), "y" = c(-.25, .1)), aes(x = x, y = y), alpha = 0.0000001) + ## add invisible points to change scales
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("Median of parameter bias across sims.")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
    xlim(c(-0.03,0.55)) + 
    #ylim(c(0,.3)) +
    scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#D95F02", "#7570B3",  "#E6AB02"),
                         labels = c("Data Augmentation", "Data Deletion-Complete","Data Deletion-Simple", "Kalman Filter", "Multiple Imputations"))
)

# parameter recovery SE ---------------------------------------------------
(gauss_paramRecovery_SE_MAR <- ggplot(data = figDat_linesMAR, aes(x = amtMiss, y = paramDiffAbsDiff_med)) +
   facet_grid(~factor(param, levels = c( "Intercept","Phi", "Beta covariates")) 
              ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC", "MNAR")), 
              scales = "free_y") + 
   geom_hline(aes(yintercept = 0), colour = "grey") + 
   #geom_errorbar(aes(ymin=paramDiff_mean - paramDiff_SD, ymax=paramDiff_mean + paramDiff_SD, color = as.factor(type)), 
   #size=0.3, width=0, position = position_dodge(width=0.03))+
   #geom_ribbon(aes(ymin = paramDiff_mean - paramDiff_SD, ymax = paramDiff_mean + paramDiff_SD, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
   geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
   #geom_point(data = data.frame("x" = c(0,.5), "y" = c(-.25, .1)), aes(x = x, y = y), alpha = 0.0000001) + ## add invisible points to change scales
   geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
   theme_classic() +
   xlab("Proportion of missing data")+ 
   theme(legend.position="top")+
   theme(legend.title=element_blank())+
   ylab("Absolute Standard Error of Parameter Recovery")+ 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
   xlim(c(-0.03,0.55)) + 
   #ylim(c(0,.3)) +
   scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#D95F02", "#7570B3",  "#E6AB02"),
                        labels = c("Data Augmentation", "Data Deletion-Complete","Data Deletion-Simple", "Kalman Filter", "Multiple Imputations"))
)
(gauss_paramRecovery_SE_MNAR <- ggplot(data = figDat_linesMNAR, aes(x = amtMiss, y = paramDiffAbsDiff_med)) +
    facet_grid(~factor(param, levels = c( "Intercept","Phi", "Beta covariates")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC", "MNAR")), 
               scales = "free_y") + 
    geom_hline(aes(yintercept = 0), colour = "grey") + 
    #geom_errorbar(aes(ymin=paramDiff_mean - paramDiff_SD, ymax=paramDiff_mean + paramDiff_SD, color = as.factor(type)), 
    #size=0.3, width=0, position = position_dodge(width=0.03))+
    #geom_ribbon(aes(ymin = paramDiff_mean - paramDiff_SD, ymax = paramDiff_mean + paramDiff_SD, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    #geom_point(data = data.frame("x" = c(0,.5), "y" = c(-.25, .1)), aes(x = x, y = y), alpha = 0.0000001) + ## add invisible points to change scales
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("Absolute Standard Error of Parameter Recovery")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
    xlim(c(-0.03,0.55)) + 
    #ylim(c(0,.3)) +
    scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#D95F02", "#7570B3",  "#E6AB02"),
                         labels = c("Data Augmentation", "Data Deletion-Complete","Data Deletion-Simple", "Kalman Filter", "Multiple Imputations"))
)

# parameter recovery coverage ---------------------------------------------
# does the confidence interval contain the true parameter? 
# calculate 95% CI for each param
figDat_cov_temp <- figDat_temp %>% 
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
  group_by(missingness, type, param, amtMiss) %>% 
  summarize(coverageNumber = sum(coverage), # the number of models that have coverage
            modelRunN = length(!is.na(coverage))# the total number of models 
            ) %>% 
  mutate(coveragePerc = coverageNumber/modelRunN)
    
figDat_covMAR <- figDat_cov %>% 
  filter(missingness %in% c("MAR: High AC", "MAR: Low AC",  "MAR: Med. AC"))
figDat_covMNAR <- figDat_cov %>% 
  filter(missingness %in% c("MNAR"))

(gauss_paramRecovery_coverage_MAR <- ggplot(data = figDat_covMAR, aes(x = amtMiss, y = coveragePerc)) +
    facet_grid(~factor(param, levels = c( "Intercept","Phi", "Beta covariates")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC", "MNAR")), 
               scales = "free_y") + 
    #geom_col(aes(x = amtMiss, y = coveragePerc, color = as.factor(type), fill = as.factor(type)), 
            # position = "dodge", alpha = .5) + 
    geom_hline(aes(yintercept = .95), colour = "grey") +
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("% of model runs where the 95% CI \n includes the simulation parameter")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
    xlim(c(-0.03,0.55)) + 
    #ylim(c(0,.3)) +
    scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#D95F02", "#7570B3",  "#E6AB02"),
                         labels = c("Data Augmentation", "Data Deletion-Complete","Data Deletion-Simple", "Kalman Filter", "Multiple Imputations")) + 
    
    scale_fill_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#D95F02", "#7570B3",  "#E6AB02"),
                         labels = c("Data Augmentation", "Data Deletion-Complete","Data Deletion-Simple", "Kalman Filter", "Multiple Imputations"))
)
(gauss_paramRecovery_coverage_MNAR <- ggplot(data = figDat_covMNAR, aes(x = amtMiss, y = coveragePerc)) +
    facet_grid(~factor(param, levels = c( "Intercept","Phi", "Beta covariates")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC", "MNAR")), 
               scales = "free_y") + 
    #geom_col(aes(x = amtMiss, y = coveragePerc, color = as.factor(type), fill = as.factor(type)), 
    # position = "dodge", alpha = .5) + 
    geom_hline(aes(yintercept = .95), colour = "grey") +
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("% of model runs where the 95% CI \n includes the simulation parameter")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
    xlim(c(-0.03,0.55)) + 
    #ylim(c(0,.3)) +
    scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#D95F02", "#7570B3",  "#E6AB02"),
                         labels = c("Data Augmentation", "Data Deletion-Complete","Data Deletion-Simple", "Kalman Filter", "Multiple Imputations")) + 
    
    scale_fill_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#D95F02", "#7570B3",  "#E6AB02"),
                        labels = c("Data Augmentation", "Data Deletion-Complete","Data Deletion-Simple", "Kalman Filter", "Multiple Imputations"))
)


# put figures together ----------------------------------------------------
(Gauss_paramRecovMAR <- ggarrange(gauss_paramRecovery_bias_MAR, gauss_paramRecovery_SE_MAR, 
                              gauss_paramRecovery_coverage_MAR, common.legend = TRUE, nrow = 1))

(Gauss_paramRecovMNAR <- ggarrange(gauss_paramRecovery_bias_MNAR, gauss_paramRecovery_SE_MNAR, 
                                  gauss_paramRecovery_coverage_MNAR, common.legend = TRUE, nrow = 1))

#(gauss_paramRecovAll <- ggarrange(Gauss_paramRecovMAR, Gauss_paramRecovMNAR, nrow = 2))

## save results
png(file = "./figures/parameterRecoveryGaussian_MAR.png", width = 12, height = 4, units = "in", res = 700)
Gauss_paramRecovMAR
dev.off()

png(file = "./figures/parameterRecoveryGaussian_MNAR.png", width = 9, height = 4, units = "in", res = 700)
Gauss_paramRecovMNAR
dev.off()
# old figures -------------------------------------------------------------

# Figure of parameter recovery (mean and sd in separate panels) -----------
# figure of means for each model type and level of missingness (with shortened x-axis)
(gauss_sim_MedsFig_trimmed <- ggplot(data = figDat_lines, aes(x = amtMiss, y = paramDiff_med)) +
   facet_grid(~factor(param, levels = c( "Intercept","Phi", "Beta covariates")) 
              ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC", "MNAR"))) + 
   geom_hline(aes(yintercept = 0), colour = "grey") + 
   #geom_errorbar(aes(ymin=paramDiff_mean - paramDiff_SD, ymax=paramDiff_mean + paramDiff_SD, color = as.factor(type)), 
   #size=0.3, width=0, position = position_dodge(width=0.03))+
   #geom_ribbon(aes(ymin = paramDiff_mean - paramDiff_SD, ymax = paramDiff_mean + paramDiff_SD, color = as.factor(type), fill = as.factor(type)), alpha = .1) +
   geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
   #geom_point(data = data.frame("x" = c(0,.5), "y" = c(-.25, .1)), aes(x = x, y = y), alpha = 0.0000001) + ## add invisible points to change scales
   geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
   theme_classic() +
   xlab("Proportion of missing data")+ 
   theme(legend.position="top")+
   theme(legend.title=element_blank())+
   ylab("Median of parameter bias across sims.")+ 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
   xlim(c(-0.03,0.55)) + 
   #ylim(c(0,.3)) +
   scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#D95F02", "#7570B3",  "#E6AB02"),
                        labels = c("Data Augmentation", "Data Deletion-Complete","Data Deletion-Simple", "Kalman Filter", "Multiple Imputations"))
)

## make a sub-dataframe that includes points for SDs that are >5 
largeSD <- figDat_lines[figDat_lines$amtMiss <= 0.5 & 
                          figDat_lines$paramDiff_SD > 5,]

# figure of SDfor each model type and level of missingness
(gauss_sim_SDFig_trimmed <- ggplot(data = figDat_lines, aes(x = amtMiss, y = paramDiff_SD)) +
    facet_grid(~factor(param, levels = c("Intercept","Phi", "Beta covariates")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC", "MNAR"))) + 
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    geom_hline(aes(yintercept = 0), colour = "grey") + theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("SD of parameter bias across sims.")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))+ 
    scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#D95F02", "#7570B3",  "#E6AB02"),
                         labels = c("Data Augmentation", "Data Deletion-Complete","Data Deletion-Simple", "Kalman Filter", "Multiple Imputations"))
  #+
  #xlim(c(-0.03,0.43)) #+ 
  #ylim(c(0,5)) #+
  #geom_point(data = largeSD, aes(x = amtMiss, y = c(3.99,3.99), color = as.factor(type)), 
  # position = position_dodge(width=0.03), pch = 8)
) 

# # figure of SD that's averaged 
#   (gauss_sim_SDFig_trimmed <- ggplot(data = figDat_lines, aes(x = amtMiss, y = SE_mean)) +
#     facet_grid(~factor(param, levels = c("Intercept","Phi", "Beta covariates")) 
#                ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC", "MNAR"))) + 
#     geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
#     geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
#     geom_hline(aes(yintercept = 0), colour = "grey") + theme_classic() +
#     xlab("Proportion of missing data")+ 
#     theme(legend.position="top")+
#     theme(legend.title=element_blank())+
#     xlim(c(0,.35)) +
#     ylim(c(0,.2)) +
#     ylab("Mean SE of parameter estimates across all sims. ")+ 
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))+ 
#     scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#D95F02", "#7570B3",  "#E6AB02"),
#                          labels = c("Data Augmentation", "Data Deletion-Complete","Data Deletion-Simple", "Kalman Filter", "Multiple Imputations"))
#   #+
#   #xlim(c(-0.03,0.43)) #+ 
#   #ylim(c(0,5)) #+
#   #geom_point(data = largeSD, aes(x = amtMiss, y = c(3.99,3.99), color = as.factor(type)), 
#   # position = position_dodge(width=0.03), pch = 8)
# ) 


# put into one figure
Gauss_paramRecov_trimmed <- ggarrange(gauss_sim_MedsFig_trimmed, gauss_sim_SDFig_trimmed, common.legend = TRUE)

## save results
png(file = "./figures/parameterRecovery_sim_Guassian_medsSD_trimmed.png", width = 9, height = 4, units = "in", res = 700)
Gauss_paramRecov_trimmed
dev.off()


#
figDat_long <- figDat_temp %>% 
  filter(param != "sigma") %>%
  filter(missingness %in% c("MAR: High AC", "MAR: Med. AC", "MAR: Low AC", "MNAR")) %>% 
  mutate(autoCor = round(autoCor, 1), 
         amtMiss = round(amtMiss, 1),
         param = replace(param, param %in% c("discharge", "light"), "Beta covariates"),
         param = replace(param, param == "intercept", "Intercept"),
         param = replace(param, param == "phi", "Phi")) %>% 
  group_by(missingness, type, param, amtMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) #%>% 
#filter(n  > 100)  

# make a figure like the one above, but without a trimmed x axis
# figure of means for each model type and level of missingness
(gauss_sim_MedsFig_reg<- ggplot(data = figDat_long, aes(x = amtMiss, y = paramDiff_med)) +
    facet_grid(~factor(param, levels = c("Intercept", "Phi", "Beta covariates")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC", "MNAR"))) + 
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
    ylab("Mean of parameter bias across sims.")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) + 
    scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#D95F02", "#7570B3",  "#E6AB02"),
                         labels = c("Data Augmentation", "Data Deletion-Complete","Data Deletion-Simple", "Kalman Filter", "Multiple Imputations"))
)

# figure of SD for each model type and level of missingness
(gauss_sim_SDFig_reg<- ggplot(data = figDat_long, aes(x = amtMiss, y = paramDiff_SD)) +
    facet_grid(~factor(param, levels = c("Intercept", "Phi", "Beta covariates")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC", "MNAR"))) + 
    geom_hline(aes(yintercept = 0), colour = "grey") + 
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    #geom_ribbon(aes(ymin = 0, ymax = paramDiff_SD, fill = as.factor(type), color = as.factor(type)), alpha = .2, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("SD of parameter bias across sims.")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))+ 
    scale_color_discrete(type = c("#1B9E77", "#66A61E", "#E7298A","#D95F02", "#7570B3",  "#E6AB02"),
                         labels = c("Data Augmentation", "Data Deletion-Complete","Data Deletion-Simple", "Kalman Filter", "Multiple Imputations"))
  
) 
# put into one figure
Gauss_paramRecov <- ggarrange(gauss_sim_MedsFig_reg, gauss_sim_SDFig_reg, common.legend = TRUE)

## save results
png(file = "./figures/parameterRecovery_sim_Guassian_medsSD.png", width = 9, height = 4, units = "in", res = 700)
Gauss_paramRecov
dev.off()
# 
# 
# 95% CI Error bar plots to show spread of complete parameter recovery data  --------
(ErrorBarPlots <- ggplot(data = figDat_long, aes(x = amtMiss, y = paramDiff_med)) +
   facet_grid(~factor(param, levels = c("Phi", "Pntercept",  "Beta covariates"))
              ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC", "MNAR"))) +
   geom_hline(aes(yintercept = 0), colour = "grey") +
   geom_errorbar(aes(ymin=paramDiff_med - 1.96*paramDiff_SD, ymax=paramDiff_med + 1.96*paramDiff_SD, color = as.factor(type)),
                 size=0.3, width=0, position = position_dodge(width=0.07))+
   geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.07)) +
   theme_classic() +
   xlab("Proportion of missing data")+
   theme(legend.position="top")+
   theme(legend.title=element_blank())+
   ylab("Mean standardized parameter estimate")+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))
)

# 
# figDat_all <- figDat_temp %>% 
#   filter(param != "sigma") %>%
#   filter(missingness %in% c("MAR: High AC", "MAR: Low AC", "MNAR")) %>% 
#   mutate(autoCor = round(autoCor, 1), 
#          amtMiss = round(amtMiss, 1)) %>% 
#   filter(amtMiss <= 0.5)
# 
# (gauss_sim_violin <- ggplot(data = figDat_temp) +
#   facet_grid(~factor(param, levels = c("intercept", "phi", "light", "discharge"))
#              ~ factor(missingness, levels = c("MAR: Low AC", "MAR_,medAutoCor", "MAR: High AC", "MNAR"))) +
#   geom_hline(aes(yintercept = 0), colour = "grey") +
#   geom_violin(aes(x = as.factor(amtMiss), y = paramDiff, color = type)) +
#   theme_classic() +
#   xlab("Proportion of missing data")+
#   theme(legend.position="top")+
#   theme(legend.title=element_blank())+
#   ylab("Mean standardized parameter estimate")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)))

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
  filter(missingness %in% c("MAR: High AC", "MAR: Low AC", "MNAR")) %>% 
  mutate(autoCor = round(autoCor, 1), 
         amtMiss = round(amtMiss, 1)) %>% 
  group_by(missingness, type, param, amtMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff),
            paramDiff_actual_mean = mean(paramDiff_actual, na.rm = TRUE),
            param_actual_sd = sd(paramDiff_actual, na.rm = TRUE)) %>% 
  filter(n  > 100)  %>% # drop combinations that have fewer than 100 observations
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

# violin plot of brms results
brmsDat <- figDat_temp %>% 
  filter(type == "brms")
ggplot(data = brmsDat) + 
  facet_grid(.~param) +
  geom_violin(aes(x = as.factor(round(amtMiss, 1)), y = paramDiff, col = round(autoCor, 1))) +
  xlab("Missingness Approach") + 
  ylab("[(true param-sim param)/abs(sim param)]") +
  ggtitle("phi parameter recovery")

