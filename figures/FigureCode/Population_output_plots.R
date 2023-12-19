#/////////////////
# This script makes figures for Population simulated data 
# 4 December 2023
#/////////////////

# Load packages
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

## read in data 
ricDat_tempA <- readRDS("./data/model_results/RickerA_resultTableAll.rds")
ricDat_tempB <- readRDS("./data/model_results/RickerB_resultTableAll.rds")
ricDat_temp <- rbind(ricDat_tempA, ricDat_tempB)
# for now, remove the simulation 176 (really tiny simulation parameters)
ricDat_temp <- ricDat_temp[ricDat_temp$SimNumber != 176,]

## for now, remove rows where simulated population went extinct (was 13% of data!)
badRows <- as.vector(sapply(ricDat_temp$drop_fits, function(x)
  ifelse(sum(str_detect(names(x), pattern = "cause"))>=1, 
             yes = "bad", 
             no = "good")
))
ricDat_temp <- ricDat_temp[rownames(ricDat_temp) %in% which(badRows == "good"),]
# rename columns 
ricDat_temp <- ricDat_temp %>% 
  rename(r_sim = r, alpha_sim = alpha, N0_sim = N0)
# put input autoCor and propMiss data into one column
ricDat_temp <- ricDat_temp %>% 
  mutate(input_args = paste0("a=", actAutoCorr, "_", "p=", actPropMiss)) %>% 
  select(-actAutoCorr, -actPropMiss)

# extract parameter information from the list columns
ricDat <- rbind(
#drop na fits
cbind(ricDat_temp[,c("SimNumber", "r", "alpha", 
                     "N0","input_args", "autoCorr", "propMiss")], 
  map_df(ricDat_temp$drop_fits, function(x) 
    ##
    if (length(names(x)) < 3) {
      data.frame(
        "type" = "dropNA",
        "r_est" =NA, 
        "alpha_est" = NA,
        "r_se" = NA,
        "alpha_se" = NA,
        "status" = "missingnessLimitReached")
    } else {
      data.frame(
        "type" = "dropNA",
        "r_est" = x$estim["r"], 
        "alpha_est" = x$estim["alpha"],
        "r_se" = x$se["r"],
        "alpha_se" = x$se["alpha"],
        "status" = "good")
    }
  ),
  data.frame("listName" = names(ricDat_temp$drop_fits))
),
#dropNA complete case fits
cbind(ricDat_temp[,c("SimNumber", "r", "alpha", 
                     "N0","input_args", "autoCorr", "propMiss")], 
      map_df(ricDat_temp$cc_fits, function(x) 
        if (length(names(x)) < 3) {
          data.frame(
            "type" = NA,
            "r_est" =NA, 
            "alpha_est" = NA,
            "r_se" = NA,
            "alpha_se" = NA,
            "status" = "missingnessLimitReached")
        } else {
          data.frame(
            "type" = "CompleteCaseDropNA",
            "r_est" = x$estim["r"], 
            "alpha_est" = x$estim["alpha"],
            "r_se" = x$se["r"],
            "alpha_se" = x$se["alpha"],
            "status" = "good")
        }
        ),
      data.frame("listName" = names(ricDat_temp$drop_fits))
),
#EM fits
cbind(ricDat_temp[,c("SimNumber", "r", "alpha", 
                     "N0","input_args", "autoCorr", "propMiss")], 
      map_df(ricDat_temp$EM_fits, function(x) 
        if (length(names(x)) < 3) {
          data.frame(
            "type" = NA,
            "r_est" =NA, 
            "alpha_est" = NA,
            "r_se" = NA,
            "alpha_se" = NA,
            "status" = "missingnessLimitReached")
        } else {
          data.frame(
            "type" = "ExpectationMaximization",
            "r_est" = x$estim["r"], 
            "alpha_est" = x$estim["alpha"],
            "r_se" = x$se["r"],
            "alpha_se" = x$se["alpha"],
            "status" = "good")
        }
      ),
      data.frame("listName" = names(ricDat_temp$drop_fits))
      
),
#DA fits
cbind(ricDat_temp[,c("SimNumber", "r", "alpha", 
                     "N0","input_args", "autoCorr", "propMiss")], 
      map_df(ricDat_temp$DA_fits, function(x) 
        if (length(names(x)) < 3) {
          data.frame(
            "type" = NA,
            "r_est" =NA, 
            "alpha_est" = NA,
            "r_se" = NA,
            "alpha_se" = NA,
            "status" = "missingnessLimitReached")
        } else {
          data.frame(
            "type" = "DataAugmentation",
            "r_est" = x$estim["r"], 
            "alpha_est" = x$estim["alpha"],
            "r_se" = x$se["r"],
            "alpha_se" = x$se["alpha"],
            "status" = "good")
        }
      ),
      data.frame("listName" = names(ricDat_temp$drop_fits))
      
),
#MI fits
cbind(ricDat_temp[,c("SimNumber", "r", "alpha", 
                     "N0","input_args", "autoCorr", "propMiss")], 
      map_df(ricDat_temp$MI_fits, function(x) 
        if (length(names(x)) < 3) {
          data.frame(
            "type" = NA,
            "r_est" =NA, 
            "alpha_est" = NA,
            "r_se" = NA,
            "alpha_se" = NA,
            "status" = "missingnessLimitReached")
        } else {
          data.frame(
            "type" = "MultipleImputations",
            "r_est" = x$estim["r"], 
            "alpha_est" = x$estim["alpha"],
            "r_se" = x$se["r"],
            "alpha_se" = x$se["alpha"],
            "status" = "good")
        }
      ),
      data.frame("listName" = names(ricDat_temp$drop_fits))
      
)
)

# remove NAs (from missingness Limit Reached issue)
ricDat <- ricDat[ricDat$status != "missingnessLimitReached",]
# get autocorrelation amount 
ricDat$actAutoCorr <- as.numeric(str_split(ricDat$listName, pattern = "_", simplify = TRUE)[,9])

#make into long data.frame 
 paramEstLong <- ricDat %>% 
  #select(SimNumber, r_sim, actAutoCorr, actPropMiss, type, r_est, r_se, status, r_paramDiff) %>% 
  pivot_longer(cols = c(r_est, alpha_est), 
               values_to = "paramEst", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
   select(-r, -alpha, -N0, -r_se, -alpha_se)
paramSimLong <- ricDat %>% 
  pivot_longer(cols = c(r, alpha), 
               values_to = "paramSim", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  select(-r_est, -alpha_est, -N0, -r_se, -alpha_se)

paramSELong <- ricDat %>% 
  pivot_longer(cols = c(r_se, alpha_se), 
               values_to = "paramSE", 
               names_to = "param", 
               names_transform = function(x) str_split(string = x, pattern = "_", simplify = TRUE)[,1]) %>% 
  select(-r, -alpha, -N0, -r_est, -alpha_est)

ricDat_long <- left_join(paramEstLong, paramSimLong) %>% 
  left_join(paramSELong)

#calculate standardized parameter estimates
ricDat_long <- ricDat_long %>%
mutate("paramDiff" = (paramEst - paramSim)/abs(paramSim))

# filter for low and high autocor
ricDat_long[ricDat_long$actAutoCorr <=0.3 & !is.na(ricDat_long$actAutoCorr), "missingness"] <- "MAR: Low AC"
ricDat_long[ricDat_long$actAutoCorr >0.3 & ricDat_long$actAutoCorr <0.6 & !is.na(ricDat_long$actAutoCorr), "missingness"] <- "MAR: Med. AC"
ricDat_long[ricDat_long$actAutoCorr  >= 0.6 & !is.na(ricDat_long$actAutoCorr), "missingness"] <- "MAR: High AC"

# save data to file for use later...
#write_rds(ricDat_long, file = "./data/model_results/ricker_sim_ModelResultsLong.rds")

# generate summary statistics for each 
ricDat_lines <- ricDat_long %>% 
  filter(missingness %in% c("MAR: High AC", "MAR: Med. AC", "MAR: Low AC")) %>% 
  mutate(autoCor = round(autoCorr, 1), 
         amtMiss = round(propMiss, 1)) %>% 
  group_by(missingness, type, param, amtMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) %>% 
  filter(n  > 50)  %>% # drop combinations that have fewer than 300 observations
  filter(amtMiss <=.5)

# Figure of parameter recovery (mean and sd in separate panels) -----------
# figure of means for each model type and level of missingness (with shortened x-axis)
(pois_sim_MedsFig_trimmed <- ggplot(data = ricDat_lines, aes(x = amtMiss, y = paramDiff_med)) +
   facet_grid(~factor(param, levels = c( "alpha", "r")) 
              ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC")),
              scales = "free_y") + 
   geom_hline(aes(yintercept = 0), colour = "grey") + 
   geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
   #geom_point(data = data.frame("x" = c(0,.5), "y" = c(-.25, .1)), aes(x = x, y = y), alpha = 0.0000001) + ## add invisible points to change scales
   geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
   #geom_line(aes(x = amtMiss, y = paramDiff_med, color = as.factor(type)), position = position_dodge(width=0.03), lty = 2) +
   theme_classic() +
   xlab("Proportion of missing data")+ 
   theme(legend.position="top", legend.title=element_blank())+
   ylab("Median of parameter bias across sims.")+ 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) + 
   scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                        labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Expectation Maximization", "Multiple Imputation"))#+
   #ylim(c(-0.85,0.55))
)
 # figure of SD for each model type and level of missingness
(pois_sim_SDFig_trimmed <- ggplot(data = ricDat_lines, aes(x = amtMiss, y = paramDiff_SD)) +
    facet_grid(~factor(param, levels = c( "alpha", "r")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC")),
               scales = "free_y") + 
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    geom_hline(aes(yintercept = 0), colour = "grey") + theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("SD of parameter bias across sims.")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) + 
    scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                         labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Expectation Maximization", "Multiple Imputation"))    
) 

 # put into one figure
pois_paramRecov_trimmed <- ggarrange(pois_sim_MedsFig_trimmed, pois_sim_SDFig_trimmed, common.legend = TRUE)

## save results
png(file = "./figures/parameterRecovery_sim_Poisson_medsSD_trimmed.png", width = 9, height = 4, units = "in", res = 700)
pois_paramRecov_trimmed
dev.off()



#
ricDat_verylong <- ricDat_long %>% 
  filter(missingness %in% c("MAR: High AC", "MAR: Med. AC", "MAR: Low AC")) %>% 
  mutate(autoCor = round(autoCorr, 1), 
         amtMiss = round(propMiss, 1)) %>% 
  group_by(missingness, type, param, amtMiss) %>% 
  summarize(paramDiff_mean = mean(paramDiff, na.rm = TRUE),
            paramDiff_med = median(paramDiff, na.rm = TRUE),
            paramDiff_SD = sd(paramDiff, na.rm = TRUE),
            n = length(paramDiff)) %>% 
  filter(n  > 100)  

# make a figure like the one above, but without a trimmed x axis
# figure of means for each model type and level of missingness
(poiss_sim_MedsFig_reg<- ggplot(data = ricDat_verylong, aes(x = amtMiss, y = paramDiff_med)) +
    facet_grid(~factor(param, levels = c("alpha", "r")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC"))) + 
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
    ylab("Median of parameter bias across sims.")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) + 
    scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                         labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Expectation Maximization", "Multiple Imputation"))
  
)
# figure of SD for each model type and level of missingness
(poiss_sim_SDFig_reg<- ggplot(data = ricDat_verylong, aes(x = amtMiss, y = paramDiff_SD)) +
    facet_grid(~factor(param, levels = c("alpha", "r")) 
               ~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC","MAR: High AC"))) + 
    geom_hline(aes(yintercept = 0), colour = "grey") + 
    geom_line(aes(color = as.factor(type)), position = position_dodge(width=0.03)) + 
    geom_point(aes(color = as.factor(type)), alpha = .8, position = position_dodge(width=0.03)) +
    #geom_ribbon(aes(ymin = 0, ymax = paramDiff_SD, fill = as.factor(type), color = as.factor(type)), alpha = .2, position = position_dodge(width=0.03)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("SD of parameter bias across sims.")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) + 
    scale_color_discrete(type = c("#66A61E","#1B9E77", "#E7298A", "#E6AB02","#7570B3"),
                         labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Expectation Maximization", "Multiple Imputation"))
  ) 
# put into one figure
poiss_paramRecov <- ggarrange(poiss_sim_MedsFig_reg, poiss_sim_SDFig_reg, common.legend = TRUE)

## save results
png(file = "./figures/parameterRecovery_sim_Poisson_medsSD.png", width = 9, height = 4, units = "in", res = 700)
poiss_paramRecov
dev.off()


ricDat_violin <- ricDat_long %>% 
  filter(missingness %in% c("MAR: High AC", "MAR: Med. AC", "MAR: Low AC")) %>% 
  mutate(autoCor = round(autoCorr, 1), 
         amtMiss = round(propMiss, 1)) %>% 
  filter(amtMiss <=0.5)


(pois_sim_violin <- ggplot(data = ricDat_violin) +
    facet_grid(~factor(param, levels = c( "alpha", "r")) 
               #~ factor(missingness, levels = c("MAR: Low AC", "MAR: Med. AC", "MAR: High AC")),
               ~ factor(type),
               scales = "free_y") + 
    geom_point(aes(x = as.factor(amtMiss), y = paramDiff, color = type), alpha = .8, position = position_dodge(width=0.07)) +
    geom_hline(aes(yintercept = 0), colour = "grey") + 
    geom_violin(aes(x = as.factor(amtMiss), y = paramDiff, color = type), draw_quantiles = c(.50)) +
    theme_classic() +
    xlab("Proportion of missing data")+ 
    theme(legend.position="top")+
    theme(legend.title=element_blank())+
    ylab("Mean standardized parameter estimate")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
    ylim(c(-8,8))
)

## save figure
png(file = "./figures/parameterRecovery_sim_Guassian_95CIs.png", width = 9, height = 4, units = "in", res = 700)
ErrorBarPlots
dev.off()
