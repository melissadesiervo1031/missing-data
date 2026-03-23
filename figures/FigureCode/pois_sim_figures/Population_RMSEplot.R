#####
# making figures for Population forecast RMSE 
#####

# load packages
library(tidyverse)
library(Metrics)
library(RColorBrewer)
library(ggpubr)
library(here)
library(janitor)
library(lazyeval)

# read in prediction data for MCAR missingness ------------------------------------------------------------
# have rmse for each run in this data.frame
ricDat_tempA <- readRDS("./data/model_results/RickerA_resultTableRev1.rds")
ricDat_tempB <- readRDS("./data/model_results/RickerB_resultTableRev1.rds")
# get MI data
ricDat_tempB_MI <- readRDS("./data/model_results/RickerB_resultTableMIRev1.rds")
ricDat_tempA_MI <- readRDS("./data/model_results/RickerA_resultTableMIRev1.rds") 
names(ricDat_tempA_MI) <- names(ricDat_tempB_MI)

# add MI to other data
ricDat_tempA <- ricDat_tempA %>% 
  cbind(ricDat_tempA_MI %>% select(forecasts_MI))

ricDat_tempB <- ricDat_tempB %>% 
  cbind(ricDat_tempB_MI  %>% select(forecasts_MI))

## get actual proportion of missingness and autocorrelation
# for version A 
missingData_A <- readRDS("./data/missingDatasets/pois_sim_randMiss_A.rds")
missingData_A_names <- lapply(X = as.vector(1:length(missingData_A)), FUN = function(x) {
  temp <- data.frame("Sim" = names(missingData_A[x]),
                     "Names" =   names(missingData_A[[x]]$y)#,
                     # "ID" = 1:length(missingData_A[[x]]$y)
  )
  return(temp)
}
) %>% 
  purrr::list_rbind()
missingData_A_names$ID <- 1:nrow(missingData_A_names)
missingData_A_names$SimNumber <- as.numeric(str_extract(missingData_A_names$Sim, pattern = "[0-9]+"))
missingData_A_names$list_AutoCorr <- str_split( missingData_A_names$Names, "_", simplify = TRUE)[,4]
missingData_A_names$list_PropMiss <- str_split( missingData_A_names$Names, "_", simplify = TRUE)[,2]
missingData_A_names[missingData_A_names$Names == "y_noMiss", c("list_AutoCorr", "list_PropMiss")] <- 0
missingData_A_names$list_AutoCorr  <- as.numeric(missingData_A_names$list_AutoCorr )
missingData_A_names$list_PropMiss <- as.numeric(missingData_A_names$list_PropMiss)

# for version B
missingData_B <- readRDS("./data/missingDatasets/pois_sim_randMiss_B.rds")
missingData_B_names <- lapply(X = as.vector(1:length(missingData_B)), FUN = function(x) {
  temp <- data.frame("Sim" = names(missingData_B[x]),
                     "Names" =   names(missingData_B[[x]]$y)#,
                     # "ID" = 1:length(missingData_B[[x]]$y)
  )
  return(temp)
}
) %>% 
  purrr::list_rbind()
missingData_B_names$ID <- 1:nrow(missingData_B_names)
missingData_B_names$SimNumber <- as.numeric(str_extract(missingData_B_names$Sim, pattern = "[0-9]+"))
missingData_B_names$list_AutoCorr <- str_split( missingData_B_names$Names, "_", simplify = TRUE)[,4]
missingData_B_names$list_PropMiss <- str_split( missingData_B_names$Names, "_", simplify = TRUE)[,2]
missingData_B_names[missingData_B_names$Names == "y_noMiss", c("list_AutoCorr", "list_PropMiss")] <- 0
missingData_B_names$list_AutoCorr  <- as.numeric(missingData_B_names$list_AutoCorr )
missingData_B_names$list_PropMiss <- as.numeric(missingData_B_names$list_PropMiss)

# add missingness info to the model results

ricDat_tempA <- ricDat_tempA %>% 
  left_join(missingData_A_names, by = c("SimNumber" = "SimNumber", "id" = "ID"))

ricDat_tempB <- ricDat_tempB %>% 
  left_join(missingData_B_names, by = c("SimNumber" = "SimNumber", "id" = "ID"))

allDat_MCAR <- rbind(ricDat_tempA, ricDat_tempB)
allDat_MCAR$missingness <- "MCAR"

# read in prediction data for MNAR missingness ------------------------------------------------------------
# have rmse for each run in this data.frame
ricDat_tempMNAR <- readRDS("./data/model_results/RickerMinMaxMissRev1.rds")

# get MI data
ricDat_MNAR_MI <- readRDS("./data/model_results/RickerMinMaxMiss_MIRev1.rds") 

# add MI to other data
ricDat_tempMNAR <- ricDat_tempMNAR %>% 
  cbind(ricDat_MNAR_MI %>% select(forecasts_MI))

## get actual proportion of missingness and autocorrelation
# for version A 
missingData_MNAR <- readRDS("./data/missingDatasets/pois_sim_minMaxMiss.rds")
missingData_MNAR_names <- lapply(X = as.vector(1:length(missingData_MNAR)), FUN = function(x) {
  temp <- data.frame("Sim" = x,
                     "Names" =   names(missingData_MNAR[[x]]$y)#,
                     # "ID" = 1:length(missingData_MNAR[[x]]$y)
  )
  return(temp)
}
) %>% 
  purrr::list_rbind()
missingData_MNAR_names$ID <- 1:nrow(missingData_MNAR_names)
missingData_MNAR_names$list_PropMiss <- str_split( missingData_MNAR_names$Names, "_", simplify = TRUE)[,2]
missingData_MNAR_names[missingData_MNAR_names$Names == "y_noMiss", c("list_PropMiss")] <- 0
missingData_MNAR_names$list_PropMiss <- as.numeric(missingData_MNAR_names$list_PropMiss)

# add missingness info to the model results

ricDat_tempMNAR <- ricDat_tempMNAR %>% 
  left_join(missingData_MNAR_names, by = c("SimNumber" = "Sim", "id" = "ID"))

allDat_MNAR <- ricDat_tempMNAR
allDat_MNAR$missingness <- "MNAR"

## add MCAR and MNAR data together
allDat_MNAR <- allDat_MNAR %>% 
  mutate(autoCorr = NA, 
         actAutoCorr = NA, 
         Sim = NA, 
         list_AutoCorr = NA) %>% 
  select(names(allDat_MCAR))
allDat <- allDat_MNAR %>% 
  rbind(allDat_MCAR) 

# prepare data for figures ------------------------------------------------
# put predictions into a 'long' format
forecasts_long <- allDat %>%
  select(-c(drop_fits, cc_fits, EM_fits, DA_fits)) %>% 
  pivot_longer(cols = c(forecast_RMSE_drop:forecasts_MI), names_to = "modelType", values_to = "RMSE") %>% 
  mutate(modelType = replace(modelType, list = c(modelType == "forecast_RMSE_drop"), values = "dropNA"),
         modelType = replace(modelType, list = c(modelType == "forecast_RMSE_cc"), values = "dropNA_cc"),
         modelType = replace(modelType, list = c(modelType == "forecast_RMSE_EM"), values = "EM"),
         modelType = replace(modelType, list = c(modelType == "forecast_RMSE_DA"), values = "DA"),
         modelType = replace(modelType, list = c(modelType == "forecasts_MI"), values = "MI"))

# if there is no missing data, remove those model runs
forecasts_long <- forecasts_long %>% 
  filter(list_PropMiss != 0) 
# actAutoCorr was incorrectly set to NA when autocor is 0, so fixing that 
# bin according to autocorrelation
forecasts_long$autocorr_binned <- NA
forecasts_long[!is.na(forecasts_long$list_AutoCorr) & forecasts_long$list_AutoCorr < 0.3 & forecasts_long$missingness == "MCAR", "autocorr_binned"] <- "low_autocorr"
forecasts_long[!is.na(forecasts_long$list_AutoCorr) & forecasts_long$list_AutoCorr >= 0.3 & forecasts_long$list_AutoCorr < 0.6 & forecasts_long$missingness == "MCAR", "autocorr_binned"] <- "med_autocorr"
forecasts_long[!is.na(forecasts_long$list_AutoCorr) & forecasts_long$list_AutoCorr >=  0.6 & forecasts_long$missingness == "MCAR", "autocorr_binned"] <- "high_autocorr"
forecasts_long$autocorr_binned <- factor(forecasts_long$autocorr_binned, 
                                         levels = c("low_autocorr",  "med_autocorr",  "high_autocorr"))
# bin according to missingness
forecasts_long$propMiss_binned <- NA
forecasts_long[!is.na(forecasts_long$list_PropMiss) & forecasts_long$list_PropMiss <= 0.3 , "propMiss_binned"] <- 0.2
forecasts_long[!is.na(forecasts_long$list_PropMiss) & forecasts_long$list_PropMiss > 0.3 & forecasts_long$list_PropMiss <= 0.5 , "propMiss_binned"] <- 0.4
forecasts_long[!is.na(forecasts_long$list_PropMiss) & forecasts_long$list_PropMiss > 0.5 , "propMiss_binned"] <- 0.6

# Plot RMSE against missingness intervals plus lines -------------------------------------------
# for this we may have to custom create segments to go with each method
group_means_MCAR <- forecasts_long %>% 
  filter(missingness == "MCAR") %>% 
  #filter(!is.infinite(RMSE)) %>% 
  mutate("RMSE_2" = RMSE,
         "RMSE_3" = RMSE) %>% 
  group_by(modelType, autocorr_binned, propMiss_binned) %>% 
  summarize(RMSE = mean(RMSE, na.rm = TRUE),
            quantile_low =  quantile(RMSE_2, na.rm = TRUE, probs = 0.25),
            quantile_high =  quantile(RMSE_3, na.rm = TRUE, probs = 0.75)) %>% 
  mutate(missingness = "MCAR") 

## for some reason there are two models w/ complete case missingness that have an RMSE of infinity? not sure why this is? (should ask Amy)
group_means_MNAR <- forecasts_long %>% 
  filter(missingness == "MNAR") %>% 
  #filter(is.finite(RMSE)) %>% 
  mutate("RMSE_2" = RMSE,
         "RMSE_3" = RMSE) %>% 
  group_by(modelType,  propMiss_binned) %>% 
  summarize(RMSE = mean(RMSE, na.rm = TRUE),
            quantile_low =  quantile(RMSE_2, na.rm = TRUE, probs = 0.25),
            quantile_high =  quantile(RMSE_3, na.rm = TRUE, probs = 0.75)) %>% 
  mutate(missingness = "MNAR", 
         autocorr_binned = NA)  %>% 
  select(names(group_means_MCAR))

## add together
RMSE_df <- group_means_MCAR %>% 
  rbind(group_means_MNAR)
# start using standard deviation for error bars
# group_SD <- tapply(forecasts_long$RMSE, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), sd, na.rm = TRUE)
# group_LQ<- group_means-1.96*group_SD
# group_HQ <- group_means+1.96*group_SD

# rmse_missingness_int <- ggplot(RMSE_df) +
#  # geom_smooth(aes(x = propMiss_binned, y = RMSE, col = modelType), method = "lm", se = FALSE) + 
#   geom_segment(aes(x=propMiss_binned,y=RMSE,ymin=quantile_low,ymax=quantile_high,col=modelType), lineend = NULL, size=0.6)+
#   geom_point(aes(x=propMiss_binned,y=RMSE,col=modelType),size=1)+
#   facet_wrap(~autocorr_binned) +
#   scale_fill_manual(values = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
#                     labels = c("Data Aug.","Data Del.-Complete", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
#   labs(
#     x = "Proportion Missingness",
#     y = "Root Mean Squared Error (RMSE)",
#     fill = "Model Type"
#   ) +
#   theme_classic() +
#   theme(
#     strip.text = element_text(size = 12),  # Customize facet labels
#     axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
#   ) #+ 
#   #ylim(c(0,25))
# 
# rmse_missingness_int

# get only medium autocorrelation
RMSE_df_med <- RMSE_df[RMSE_df$autocorr_binned=="med_autocorr" | RMSE_df$missingness == "MNAR",]

(rmse_NoLineErrorBar <- ggplot(data = RMSE_df_med) +
    facet_grid(.~missingness) +
    geom_linerange(aes(x = propMiss_binned, ymin = quantile_low, ymax = quantile_high, color = modelType), alpha = 1, position = position_dodge(width = .1)) +
    geom_point(aes(x = propMiss_binned, y = RMSE, color = modelType), alpha = 1, position = position_dodge(width = .1)) +
    #geom_smooth(aes(x = propMiss_bin, y = RMSE_mean, col = type), method = "lm", se = FALSE) +
    theme_classic() +
    ylab("Root Mean Square Error (RMSE)") +
    xlab("Proportion of Missing Data") + 
    #ylim(c(0,1.25)) + 
    scale_x_continuous(breaks=c(0.2,0.4, 0.6)) +
    
    scale_color_discrete(type = c("#CC79A7", "#E69F00", "#D55E00", "#8c8c8c", "#009E73"),
                         labels = c("Data Augmentation", "Data Deletion-Simple", "Data Deletion-Complete", "Expectation Maximization", "Multiple Imputation")
    )  +
    guides(col = guide_legend(title = "Model Type", position = "top", direction = "vertical", nrow = 2))
)

png(file = "./figures/RMSE_FullFigure_NoLineWithErrorBar_poisson_SIM.png", width = 6, height = 6, units = "in", res = 700)
rmse_NoLineErrorBar
dev.off()
# 
# # right factor order
# RMSE_df$modelType=factor(RMSE_df$modelType,levels=c("dropNA","dropCC","MI","EM","DA"))
# custom_seg$modelType=factor(custom_seg$modelType,levels=c("dropNA","dropCC","MI","EM","DA"))
# 
# rmse_missingness_int_med <- ggplot(RMSE_df_med) +
#   #geom_smooth(aes(x = propMiss, y = RMSE, col = modelType), method = "lm", se = FALSE) + 
#   geom_segment(data=custom_seg,aes(x=x,y=y,xend=xend,yend=yend,col=modelType),size=0.6)+
#   geom_point(data=custom_seg,aes(x=x,y=means1,col=modelType),size=1)+
#   scale_x_continuous(breaks = c(0.2,0.4,0.6)) +
#   scale_color_discrete(type = c("#E69F00", "#D55E00","#009E73", "#BBBBBB","#CC79A7"),
#                        labels = c("Data Deletion-Simple", "Data Deletion-Complete", "Multiple Imputation", "Expectation Maximization", "Data Augmentation")) +
#   labs(
#     x = "Proportion Missing",
#     y = "Root Mean Square Error (RMSE)",
#     color = "Model Type"
#   ) +
#   theme_classic() +
#   theme(legend.position="top")+
#   theme(legend.title=element_blank())+
#   guides(color = guide_legend(nrow = 2))+
#   theme(
#     strip.text = element_text(size = 12),  # Customize facet labels
#     axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
#   )
# 
# rmse_missingness_int_med
# 
# # add the wytham tits to the plot
# wytham <- read.csv(here("data/Wytham_tits.csv"), header=TRUE)
# 
# wytham_tits = ggplot() + 
#   geom_rect(aes(xmin = 2008.5, xmax = 2018.5, 
#                 ymin = 100, ymax = 500), fill = "grey80") +
#   geom_line(data = wytham, aes(x = Year, y = Broods)) + 
#   geom_point(data = wytham, aes(x = Year, y = Broods),fill="gray",color="black",size=3,shape=21) + 
#   theme_classic() + 
#   labs(x = "Year", y = "Number of Broods")
# wytham_tits
# 
# png(file = "./figures/RMSE_pois_combined.png", width = 6, height = 7, units = "in", res = 700)
# ggarrange(wytham_tits,rmse_missingness_int_med,
#           labels = c("A", "B"),
#           ncol=1,nrow=2,
#           heights = c(1, 2))
# dev.off()
# 
# # Figure of RMSE variation  -----------------------------------------------
# # calculate the width of the IQR for each point shown in the previous figure
# RMSE_errorBar <- custom_seg %>% 
#   mutate(IQR_width = yend - y)
# 
# (rmseWidth_fig <- ggplot(data = RMSE_errorBar) +
#     geom_linerange(aes(x = xend, ymin = 0, ymax = IQR_width, color = modelType), alpha = 1, position = position_dodge(width = .1), lwd = 2) +
#     theme_classic() +
#     ylab("Width of RMSE Inter-Quartile Range") +
#     xlab("Proportion of Missing Data") + 
#     #ylim(c(0,1.25)) + 
#     scale_x_continuous(breaks=c(0.2,0.4, 0.6)) +
#     scale_color_discrete(type = c("#E69F00", "#D55E00","#009E73", "#BBBBBB","#CC79A7"),
#                          labels = c("Data Deletion-Simple", "Data Deletion-Complete", "Multiple Imputation", "Expectation Maximization", "Data Augmentation")) +
#     guides(col = guide_legend(title = "Model Type", position = "top", direction = "vertical", nrow = 3))
# )
# 
# 
# png(file = "./figures/RMSE_IQR_width_poisson_auSable.png", width =5, height = 5, units = "in", res = 700)
# rmseWidth_fig
# dev.off()
# 
# 
# # Make figure of CI coverage -------------------------------------------
# # deprecated for now...
# # Error/confidence intervals should expand over time like for simulations with temporal autocorrelation
# # This is because we are basing each next data point prediction on the last data point prediction
# # We will therefore need to re-calculate the CI- I think we would do it in the modelRuns_rickerForecasts.R script
# # Use the lower and upper estimates to expand a lower and upper CI over time
# # We didn't need to do this with the simulated data because it was just measuring coverage of a simple CI
# # Will need to add extra prediction columns for each
# 
# # make coverage into proportion
# # RMSE_df$coverage=RMSE_df$coverage/10
# # 
# # # simple scatterplot- kind of ugly honestly ----------------------------------------------
# # coverage_plot <- ggplot(RMSE_df) +
# #   geom_point(aes(x = propMiss, y = coverage, col = modelType), alpha = .5) +
# #   geom_smooth(aes(x = propMiss, y = coverage, col = modelType), method = "lm", se = FALSE) + 
# #   facet_wrap(~autocorr_binned, ) +
# #   scale_color_discrete(type = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
# #                        labels = c("Data Aug.","Data Del.-Complete",  "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
# #   theme_classic()
# # 
# # coverage_plot
# 
# # simple boxplot- even more ugly honestly ----------------------------------------------
# # coverage_box <- ggplot(RMSE_df) +
# #   geom_boxplot(aes(x = factor(propMissCat), y = coverage, fill = modelType), alpha = 0.7) +
# #   facet_wrap(~autocorr_binned) +
# #   scale_fill_manual(values = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
# #                     labels = c("Data Aug.","Data Del.-Complete", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
# #   labs(
# #     x = "Proportion Missingness",
# #     y = "Coverage",
# #     fill = "Model Type"
# #   ) +
# #   theme_classic() +
# #   theme(
# #     strip.text = element_text(size = 12),  # Customize facet labels
# #     axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
# #   )
# # 
# # coverage_box
# 
# # Plot coverage against missingness intervals plus lines -------------------------------------------
# # for this we may have to custom create segments to go with each method
# # group_means <- tapply(RMSE_df$coverage, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), mean, na.rm = TRUE)
# # group_LQ<- tapply(RMSE_df$coverage, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), quantile, na.rm = TRUE, probs=0.25)
# # group_HQ <- tapply(RMSE_df$coverage, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), quantile, na.rm = TRUE, probs=0.75)
# # custom_seg=data.frame(
# #   x=c(rep(as.numeric(colnames(group_LQ[,1,])),each=nrow(group_LQ[,1,])),rep(as.numeric(colnames(group_LQ[,2,])),each=nrow(group_LQ[,2,])),rep(as.numeric(colnames(group_LQ[,3,])),each=nrow(group_LQ[,3,])))+rep(c(-0.024,-0.012,0,0.012,0.024),times=12),
# #   xend=c(rep(as.numeric(colnames(group_LQ[,1,])),each=nrow(group_LQ[,1,])),rep(as.numeric(colnames(group_LQ[,2,])),each=nrow(group_LQ[,2,])),rep(as.numeric(colnames(group_LQ[,3,])),each=nrow(group_LQ[,3,])))+rep(c(-0.024,-0.012,0,0.012,0.024),times=12),
# #   y=c(as.vector(group_LQ[,1,]),as.vector(group_LQ[,2,]),as.vector(group_LQ[,3,])),
# #   yend=c(as.vector(group_HQ[,1,]),as.vector(group_HQ[,2,]),as.vector(group_HQ[,3,])),
# #   autocorr_binned=rep(c("low_autocorr","med_autocorr","high_autocorr"),each=nrow(group_LQ[,1,])*ncol(group_LQ[,1,])),
# #   modelType=rep(rownames(group_LQ[,1,]),12),
# #   means1=c(as.vector(group_means[,1,]),as.vector(group_means[,2,]),as.vector(group_means[,3,]))
# # )
# # 
# # 
# # RMSE_df$autocorr_binned <- factor(RMSE_df$autocorr_binned, levels = c("low_autocorr","med_autocorr","high_autocorr"))
# # custom_seg$autocorr_binned <- factor(custom_seg$autocorr_binned, levels = c("low_autocorr","med_autocorr","high_autocorr"))
# # 
# # 
# # coverage_int <- ggplot(RMSE_df) +
# #   geom_smooth(aes(x = propMiss, y = coverage, col = modelType), method = "lm", se = FALSE) + 
# #   geom_segment(data=custom_seg,aes(x=x,y=y,xend=xend,yend=yend,col=modelType),size=0.6)+
# #   geom_point(data=custom_seg,aes(x=x,y=means1,col=modelType),size=1)+
# #   facet_wrap(~autocorr_binned) +
# #   scale_fill_manual(values = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
# #                     labels = c("Data Aug.","Data Del.-Complete", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
# #   labs(
# #     x = "Proportion Missingness",
# #     y = "Coverage",
# #     fill = "Model Type"
# #   ) +
# #   theme_classic() +
# #   theme(
# #     strip.text = element_text(size = 12),  # Customize facet labels
# #     axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
# #   )
# # 
# # coverage_int
# 
# 
# # Make figure of ALTERNATE CI coverage -------------------------------------------
# # deprecated for now
# # Error/confidence intervals should expand over time like for simulations with temporal autocorrelation
# # This is because we are basing each next data point prediction on the last data point prediction
# # We will therefore need to re-calculate the CI- I think we would do it in the modelRuns_rickerForecasts.R script
# # Use the lower and upper estimates to expand a lower and upper CI over time
# # We didn't need to do this with the simulated data because it was just measuring coverage of a simple CI
# # Will need to add extra prediction columns for each
# 
# # make coverage into proportion
# # RMSE_df$coverage_opt2=RMSE_df$coverage_opt2/10
# # 
# # # simple scatterplot- kind of ugly honestly ----------------------------------------------
# # coverage_opt2_plot <- ggplot(RMSE_df) +
# #   geom_point(aes(x = propMiss, y = coverage_opt2, col = modelType), alpha = .5) +
# #   geom_smooth(aes(x = propMiss, y = coverage_opt2, col = modelType), method = "lm", se = FALSE) + 
# #   facet_wrap(~autocorr_binned, ) +
# #   scale_color_discrete(type = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
# #                        labels = c("Data Aug.","Data Del.-Complete",  "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
# #   theme_classic()
# # 
# # coverage_opt2_plot
# # 
# # # simple boxplot- even more ugly honestly ----------------------------------------------
# # coverage_opt2_box <- ggplot(RMSE_df) +
# #   geom_boxplot(aes(x = factor(propMissCat), y = coverage_opt2, fill = modelType), alpha = 0.7) +
# #   facet_wrap(~autocorr_binned) +
# #   scale_fill_manual(values = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
# #                     labels = c("Data Aug.","Data Del.-Complete", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
# #   labs(
# #     x = "Proportion Missingness",
# #     y = "Coverage",
# #     fill = "Model Type"
# #   ) +
# #   theme_classic() +
# #   theme(
# #     strip.text = element_text(size = 12),  # Customize facet labels
# #     axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
# #   )
# # 
# # coverage_opt2_box
# # 
# # # Plot coverage_opt2 against missingness intervals plus lines -------------------------------------------
# # # for this we may have to custom create segments to go with each method
# # group_means <- tapply(RMSE_df$coverage_opt2, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), mean, na.rm = TRUE)
# # group_LQ<- tapply(RMSE_df$coverage_opt2, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), quantile, na.rm = TRUE, probs=0.25)
# # group_HQ <- tapply(RMSE_df$coverage_opt2, list(RMSE_df$modelType,RMSE_df$autocorr_binned,RMSE_df$propMissCat), quantile, na.rm = TRUE, probs=0.75)
# # custom_seg=data.frame(
# #   x=c(rep(as.numeric(colnames(group_LQ[,1,])),each=nrow(group_LQ[,1,])),rep(as.numeric(colnames(group_LQ[,2,])),each=nrow(group_LQ[,2,])),rep(as.numeric(colnames(group_LQ[,3,])),each=nrow(group_LQ[,3,])))+rep(c(-0.024,-0.012,0,0.012,0.024),times=12),
# #   xend=c(rep(as.numeric(colnames(group_LQ[,1,])),each=nrow(group_LQ[,1,])),rep(as.numeric(colnames(group_LQ[,2,])),each=nrow(group_LQ[,2,])),rep(as.numeric(colnames(group_LQ[,3,])),each=nrow(group_LQ[,3,])))+rep(c(-0.024,-0.012,0,0.012,0.024),times=12),
# #   y=c(as.vector(group_LQ[,1,]),as.vector(group_LQ[,2,]),as.vector(group_LQ[,3,])),
# #   yend=c(as.vector(group_HQ[,1,]),as.vector(group_HQ[,2,]),as.vector(group_HQ[,3,])),
# #   autocorr_binned=rep(c("low_autocorr","med_autocorr","high_autocorr"),each=nrow(group_LQ[,1,])*ncol(group_LQ[,1,])),
# #   modelType=rep(rownames(group_LQ[,1,]),12),
# #   means1=c(as.vector(group_means[,1,]),as.vector(group_means[,2,]),as.vector(group_means[,3,]))
# # )
# # 
# # 
# # RMSE_df$autocorr_binned <- factor(RMSE_df$autocorr_binned, levels = c("low_autocorr","med_autocorr","high_autocorr"))
# # custom_seg$autocorr_binned <- factor(custom_seg$autocorr_binned, levels = c("low_autocorr","med_autocorr","high_autocorr"))
# # 
# # 
# # coverage_opt2_int <- ggplot(RMSE_df) +
# #   geom_smooth(aes(x = propMiss, y = coverage_opt2, col = modelType), method = "lm", se = FALSE) + 
# #   geom_segment(data=custom_seg,aes(x=x,y=y,xend=xend,yend=yend,col=modelType),size=0.6)+
# #   geom_point(data=custom_seg,aes(x=x,y=means1,col=modelType),size=1)+
# #   facet_wrap(~autocorr_binned) +
# #   scale_fill_manual(values = c("#1B9E77","#66A61E", "#E7298A", "#E6AB02","#7570B3"),
# #                     labels = c("Data Aug.","Data Del.-Complete", "Data Del.-Simple", "Expectation Max.", "Multiple Imp.")) +
# #   labs(
# #     x = "Proportion Missingness",
# #     y = "Coverage",
# #     fill = "Model Type"
# #   ) +
# #   theme_classic() +
# #   theme(
# #     strip.text = element_text(size = 12),  # Customize facet labels
# #     axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
# #   )
# # 
# # coverage_opt2_int
# # 
# # 
# # 
# # # Make figure of predictions ----------------------------------------------
# # forecasts_long$propMiss_binned <- NA
# # forecasts_long[forecasts_long$propMiss <0.2, "propMiss_binned"] <- "0-0.2"
# # forecasts_long[c(forecasts_long$propMiss >=0.2 & forecasts_long$propMiss < 0.4), "propMiss_binned"] <- "0.2-0.4"
# # forecasts_long[forecasts_long$propMiss >=0.4 & forecasts_long$propMiss < 0.6, "propMiss_binned"] <- "0.4-0.6"
# # forecasts_long[forecasts_long$propMiss >=0.6 & forecasts_long$propMiss < 0.8, "propMiss_binned"] <- "0.6-0.8"
# # forecasts_long[forecasts_long$propMiss >=0.8, "propMiss_binned"] <- "0.8-1"
# # 
# # forecasts_long$propMissCat=NA
# # forecasts_long$propMissCat[which(forecasts_long$propMiss==0)]=0
# # forecasts_long$propMissCat[which(forecasts_long$propMiss>=0.15&forecasts_long$propMiss<=0.25)]=0.2
# # forecasts_long$propMissCat[which(forecasts_long$propMiss>=0.35&forecasts_long$propMiss<=0.45)]=0.4
# # forecasts_long$propMissCat[which(forecasts_long$propMiss>=0.55&forecasts_long$propMiss<=0.65)]=0.6
# # 
# # forecasts_avg <- forecasts_long %>% 
# #   #filter(Year > 2013) %>% 
# #   group_by(Year, missingness, autocorr_binned, propMissCat) %>% 
# #   summarize(Estimate_mean = mean(Estimate),
# #             Estimate_sd = sd(Estimate)) #%>% 
# # #,Est.Error_mean = mean(Est.Error)) %>% 
# # #rename(propMiss = propMiss_binned)
# # 
# # # calculate low, med, and high autocorr
# # ### set color palettes ------------------------------------------------------
# # palCont <- rev(RColorBrewer::brewer.pal(n = 4, name = "Reds"))
# # 
# # (forecastFig_p <- ggplot() + 
# #     facet_grid(.~as.factor(missingness) ~ as.factor(autocorr_binned)) +
# #     geom_line(data = realDat, aes(x = Year, y = Broods)) + 
# #     geom_line(data = forecasts_avg, aes(x = Year, y = Estimate_mean, col = as.factor(propMissCat), group = propMissCat), alpha = .8) +
# #     theme_classic() +
# #     scale_color_viridis_d(option = "plasma", name = "prop. missing") +
# #     xlim(2005,2019)
# # )
# # #
# # # save figure
# # png(file = "./figures/forecastAccuracy_poisson.png", width = 9, height = 6, units = "in", res = 700)
# # forecastFig_p
# # dev.off()
# # 
# # 
# # # figure of mean predictions across all amounts missingness for each type  --------
# # forecasts_fig <- forecasts_avg
# # #forecasts_fig <- forecasts_avg %>% 
# # #  filter(propMiss <=0.5)
# # 
# # ggplot() + 
# #   facet_grid(.~as.factor(autocorr_binned)~ as.factor(propMissCat)) +
# #   #geom_ribbon(data = second_fig, aes(x = date, ymin = Estimate_mean - 1.96 * Est.Error_mean, ymax = Estimate_mean + 1.96 * Est.Error_mean, fill = type, group =type), alpha = .3) +
# #   geom_line(data = realDat, aes(x = Year, y = Broods)) + 
# #   geom_line(data = forecasts_fig, aes(x = Year, y = Estimate_mean, col = missingness, group = missingness)) +
# #   theme_classic() +
# #   xlim(2009,2019)
# # 
# # # figure of mean predictions across all amounts of missingness for each type  --------
# # forecasts_fig3 <- forecasts_long %>% 
# #   group_by(Year, missingness, autocorr_binned) %>% 
# #   summarize(Estimate_mean = mean(Estimate),
# #             Estimate_sd = sd(Estimate)) 
# # 
# # ggplot() + 
# #   facet_grid(.~as.factor(missingness)~ as.factor(autocorr_binned)) +
# #   geom_ribbon(data = forecasts_fig3, aes(x = Year, ymin = Estimate_mean - 1.96 * Estimate_sd, ymax = Estimate_mean + 1.96 * Estimate_sd, 
# #                                          fill = missingness, group = missingness), alpha = .3) +
# #   geom_line(data = realDat, aes(x = Year, y = Broods)) + 
# #   geom_line(data = forecasts_fig3, aes(x = Year, y= Estimate_mean, group = missingness), alpha = .3) +
# #   theme_classic() +
# #   xlim(2009,2019)
# 
# 
#/////////