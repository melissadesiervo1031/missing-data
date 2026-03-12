#/////////////////
# This script makes figures for GPP simulated data 
# 30 October 2023
#/////////////////

# Load packages
library(tidyverse)
library(ggpubr)
library(forecast)
library(RColorBrewer)
library(ggh4x)

########
figDat_temp <- readRDS("./data/model_results/gauss_real_ModelResults.rds") #

#figDat_temp <- figDat_temp[!(figDat_temp$simName %in% c(376, 831, 816, 461, 808, 129, 366, 208, 385)),]

# filter for low and high autocor\
figDat_temp[figDat_temp$missingness=="MAR" & figDat_temp$autoCor <=0.3, "missingness"] <- "MCAR: Low AC"
figDat_temp[figDat_temp$missingness=="MAR" & figDat_temp$autoCor > 0.3 & figDat_temp$autoCor < 0.6, "missingness"] <- "MCAR: Med. AC"
figDat_temp[figDat_temp$missingness=="MAR" & figDat_temp$autoCor >=  0.6, "missingness"] <- "MCAR: High AC"

# reorder factor levels for plotting ##
figDat_temp <- figDat_temp %>% 
  mutate(type=fct_relevel(type,c("Data Deletion Simple", "Data Deletion CC" ,"Kalman filter", "Multiple imputations","brms")))#%>% 
 # mutate(missingness = replace(missingness, missingprop_autocor == "y_noMiss", "MCAR: Med. AC"))


# figDat_temp %>%
#   filter(missingprop_autocor != "y_noMiss",
#          parameters != "sigma") %>%
#   filter(parameters == "phi") %>%
#   ggplot() +
#   facet_wrap(.~missingness) +
#   geom_point(aes(x = amtMiss, y = paramDiff, col = type), alpha = .2) +
#   geom_smooth(aes(x = amtMiss, y = paramDiff, col = type)) %>% 
  
figDat_lines <- figDat_temp %>% 
  filter(parameters != "sigma",
         missingprop_autocor != "y_noMiss"
  ) %>%
  # try removing simulations that have a parameterseter that is very, very small
  # filter( parameters_simVal > 0.05) %>% 
  mutate(autoCor = round(autoCor, 1), 
         amtMiss = round(amtMiss,1),
         amtMiss = replace(amtMiss, amtMiss <=0.3, 0.2),
         amtMiss = replace(amtMiss, amtMiss > 0.3 & amtMiss <=0.5, 0.4),
         amtMiss = replace(amtMiss, amtMiss > 0.5, 0.6),
         parameters = replace(parameters, parameters %in% c("discharge", "light"), "Beta covariates"),
         parameters = replace(parameters, parameters == "intercept", "Intercept"),
         parameters = replace(parameters, parameters == "phi", "Phi")
  ) %>% 
  group_by(missingness, type, parameters, amtMiss) %>% 
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
  
  #filter(n  > 50)  %>% # drop combinations that have fewer than 300 observations
  filter(amtMiss <=.65, 
         parameters != "Intercept")



figDat_lines2 <- figDat_lines %>% filter(missingness=="MCAR: Med. AC"|missingness=="MNAR")  %>%  filter(parameters=="Beta covariates"| parameters=="Phi")

figDat_lines2 %>% 
  table()
#figDat_temp %>%
 # filter(missingness %in% c("MCAR: Med. AC", "MNAR")) %>% 
  #filter(parameters %in% c("Intercept", "intercept")) %>%
  # ggplot() +
  # facet_grid(rows = vars(parameters),
  #   cols = vars(missingness), scales = "free") +
  # geom_point(aes(x = amtMiss, y = param_value, col = as.factor(type)), alpha = .5) + 
  # geom_smooth(aes(x = amtMiss, y = param_value, col = as.factor(type)))
# parameter recovery bias -------------------------------------------------
(gauss_paramRecovery_bias <- ggplot(data = figDat_lines2, aes(x = amtMiss, y = paramDiff_med)) +
   ggh4x::facet_grid2(~factor(parameters, levels = c("Phi", "Beta covariates"))
                      ~ factor(missingness, levels = c("MCAR: Med. AC", "MNAR"), labels =c("Missing Completely at Random", "Missing NOT at Random")), scales = "free_y")+
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
   ylab("Median Error")+ 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))+
   #xlim(c(0.15,0.65)) + 
   #ylim(c(0,.3)) +
   scale_color_discrete(type = c( "#CC79A7",  "#D55E00","#E69F00","#0072B2", "#009E73"),
                        labels = c( "Data Augmentation",  "Data Deletion-Complete","Data Deletion-Simple","Kalman Filter","Multiple Imputation"))
 # scale_color_discrete(type = c("#E69F00", "#D55E00","#009E73","#0072B2", "#CC79A7"),
 #                      labels = c("Data Deletion-Simple", "Data Deletion-Complete","Multiple Imputation","Kalman Filter", "Data Augmentation"))
)



# parameter recovery SE ---------------------------------------------------
(gauss_paramRecovery_SE <- ggplot(data = figDat_lines2, aes(x = amtMiss, y = paramDiffAbsDiff_med)) +
   ggh4x::facet_grid2(~factor(parameters, levels = c("Phi", "Beta covariates")) 
                      ~ factor(missingness, levels = c("MCAR: Med. AC", "MNAR"), labels =c("Missing Completely at Random", "Missing NOT at Random")), scales = "free_y")+
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
   ylab("Median Absolute Error")+ 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
   #xlim(c(0.15,0.65))+ 
   #ylim(c(0,.3)) +
   scale_color_discrete(type = c( "#CC79A7",  "#D55E00","#E69F00","#0072B2", "#009E73"),
                        labels = c( "Data Augmentation",  "Data Deletion-Complete","Data Deletion-Simple","Kalman Filter","Multiple Imputation"))
 
)


# parameter recovery coverage ---------------------------------------------
# does the confidence interval contain the true parameter? 
# calculate 95% CI for each param
figDat_cov_temp <- figDat_temp %>% 
  mutate(CI95_lower = param_value - 1.96*param_se,
         CI95_upper = param_value + 1.96*param_se) %>% 
  filter(parameters != "sigma") %>% 
  filter(!is.na(param_se))# randomly there are some model runs that don't have SE? 

# is the true parameter within the 95% CI? 
figDat_cov_temp$coverage <- c(figDat_cov_temp$param_simVal >= figDat_cov_temp$CI95_lower & 
                                figDat_cov_temp$param_simVal <= figDat_cov_temp$CI95_upper)

## count the # of models w/ and without coverage for each bin of missingness and autocorrelation
figDat_cov <- figDat_cov_temp %>% 
  filter(parameters != "sigma",
         parameters != "intercept",
         amtMiss <=.65,) %>% 
  mutate(autoCor = round(autoCor, 1), 
         #amtMiss = round(amtMiss, 1),
         amtMiss = replace(amtMiss, amtMiss <=0.3 & amtMiss > 0, 0.2),
         amtMiss = replace(amtMiss, amtMiss > 0.3 & amtMiss <=0.5, 0.4),
         amtMiss = replace(amtMiss, amtMiss > 0.5, 0.6),
         parameters = replace(parameters, parameters %in% c("discharge", "light"), "Beta covariates"),
         #parameters = replace(parameters, parameters == "intercept", "Intercept"),
         parameters = replace(parameters, parameters == "phi", "Phi")
  ) %>% 
  group_by(missingness, type, parameters, amtMiss) %>% 
  summarize(coverageNumber = sum(coverage), # the number of models that have coverage
            modelRunN = length(!is.na(coverage))# the total number of models 
  ) %>% 
  mutate(coveragePerc = coverageNumber/modelRunN)


figDat_cov2 <- figDat_cov %>% filter(missingness=="MCAR: Med. AC"|missingness=="MNAR")  %>%  filter(parameters=="Beta covariates"| parameters=="Phi")



(gauss_paramRecovery_coverage <- ggplot(data = figDat_cov2, aes(x = amtMiss, y = coveragePerc)) +
    
    #geom_col(aes(x = amtMiss, y = coveragePerc, color = as.factor(type), fill = as.factor(type)), 
    # position = "dodge", alpha = .5) + 
    ggh4x::facet_grid2(~factor(parameters, levels = c("Phi", "Beta covariates")) 
                       ~ factor(missingness, levels = c("MCAR: Med. AC", "MNAR"), labels =c("Missing Completely at Random", "Missing NOT at Random")), scales = "free_y")+
    geom_hline(aes(yintercept = .95), colour = "grey", linetype="dotted") +
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("Coverage")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
    #xlim(c(0.15,0.65)) + 
    #ylim(c(0,.3)) +
    scale_color_discrete(type = c( "#CC79A7",  "#D55E00","#E69F00","#0072B2", "#009E73"),
                         labels = c( "Data Augmentation",  "Data Deletion-Complete","Data Deletion-Simple","Kalman Filter","Multiple Imputation"))
  
)


# put figures together long version ----------------------------------------------------


gauss_paramRecovery_bias2<-gauss_paramRecovery_bias+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), legend.text = element_text(size=7))

gauss_paramRecovery_SE2<-gauss_paramRecovery_SE+ theme(strip.text.x = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),legend.text = element_text(size=7))

gauss_paramRecovery_coverage2<-gauss_paramRecovery_coverage+ theme(strip.text.x = element_blank(),legend.text = element_text(size=7))

(Gauss_paramRecov_MCARMNARlong <- ggarrange(gauss_paramRecovery_bias2, gauss_paramRecovery_SE2, 
                                            gauss_paramRecovery_coverage2, common.legend = TRUE, nrow = 3, align = "v"))

#(gauss_paramRecovAll <- ggarrange(Gauss_paramRecovMAR, Gauss_paramRecovMNAR, nrow = 2))

## save results
png(file = "./figures/parameterRecoveryGaussian_REAL_MARMNARlong.png", width = 6.5, height = 8, units = "in", res = 700)
Gauss_paramRecov_MCARMNARlong
dev.off()


# put figures together wide version ----------------------------------------------------

gauss_paramRecovery_bias<-gauss_paramRecovery_bias+ theme(strip.text.y = element_blank())

gauss_paramRecovery_SE<-gauss_paramRecovery_SE+ theme(strip.text.y = element_blank())

(Gauss_paramRecov_MARMNAR <- ggarrange(gauss_paramRecovery_bias, gauss_paramRecovery_SE, 
                                       gauss_paramRecovery_coverage, common.legend = TRUE, nrow = 1))

#(gauss_paramRecovAll <- ggarrange(Gauss_paramRecovMAR, Gauss_paramRecovMNAR, nrow = 2))

## save results
png(file = "./figures/parameterRecoveryGaussian_REAL_MARMNAR.png", width = 13.5, height = 4, units = "in", res = 700)
Gauss_paramRecov_MARMNAR
dev.off()



