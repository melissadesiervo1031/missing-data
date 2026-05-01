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
forecasts_long[!is.na(forecasts_long$list_AutoCorr) & forecasts_long$list_AutoCorr  <= 0.25  & forecasts_long$missingness == "MCAR", "autocorr_binned"] <- "low_autocorr"
forecasts_long[!is.na(forecasts_long$list_AutoCorr) & forecasts_long$list_AutoCorr >0.25 & forecasts_long$list_AutoCorr < 0.65 & forecasts_long$missingness == "MCAR", "autocorr_binned"] <- "med_autocorr"
forecasts_long[!is.na(forecasts_long$list_AutoCorr) & forecasts_long$list_AutoCorr >=  0.65 & forecasts_long$missingness == "MCAR", "autocorr_binned"] <- "high_autocorr"
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
            RMSE_med = median(RMSE_2, na.rm = TRUE),
            quantile_low =  quantile(RMSE_2, na.rm = TRUE, probs = 0.25),
            quantile_high =  quantile(RMSE_3, na.rm = TRUE, probs = 0.75),
            RMSE_n = sum(!is.na(RMSE_2)))%>% 
  mutate(missingness = "MCAR") 

## for some reason there are two models w/ complete case missingness that have an RMSE of infinity? not sure why this is? (should ask Amy)
group_means_MNAR <- forecasts_long %>% 
  filter(missingness == "MNAR") %>% 
  #filter(is.finite(RMSE)) %>% 
  mutate("RMSE_2" = RMSE,
         "RMSE_3" = RMSE) %>% 
  group_by(modelType,  propMiss_binned) %>% 
  summarize(RMSE = mean(RMSE, na.rm = TRUE),
            RMSE_med = median(RMSE_2, na.rm = TRUE),
            quantile_low =  quantile(RMSE_2, na.rm = TRUE, probs = 0.25),
            quantile_high =  quantile(RMSE_3, na.rm = TRUE, probs = 0.75),
            RMSE_n = sum(!is.na(RMSE_2))) %>% 
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
    ggh4x::facet_grid2(
                       ~ factor(missingness, levels = c("MCAR", "MNAR"), labels =c("'Missing Completely at Random'", "'Missing Not at Random'")),
                       labeller =  label_parsed) +
    geom_linerange(aes(x = propMiss_binned, ymin = quantile_low, ymax = quantile_high, color = modelType), alpha = 1, position = position_dodge(width = .1)) +
    geom_point(aes(x = propMiss_binned, y = RMSE_med, color = modelType), alpha = 1, position = position_dodge(width = .1)) +
    #geom_text(aes(x = propMiss_binned, y = 0.8, label = paste0("N=",RMSE_n), group = modelType),  hjust = 0.2, angle =90, position = position_dodge(width = .1) ) + 
    #geom_smooth(aes(x = propMiss_bin, y = RMSE_mean, col = type), method = "lm", se = FALSE) +
    theme_classic() +
    ylab("Root Mean Square Error (RMSE)") +
    xlab("Proportion of Missing Data") + 
    #ylim(c(0,1.25)) + 
    scale_x_continuous(breaks=c(0.2,0.4, 0.6)) +
    
    scale_color_discrete(type = c("#CC79A7", "#E69F00", "#D55E00", "#8c8c8c", "#009E73"),
                         labels = c("Data Augmentation", "Data Deletion-Simple", "Data Deletion-Complete", "Expectation Maximization", "Multiple Imputation")
    )  +
    guides(col = guide_legend(title = "Model Type", position = "right", direction = "vertical", nrow = 5))
)

png(file = "./figures/RMSE_FullFigure_NoLineWithErrorBar_poisson_SIM.png", width = 6, height = 6, units = "in", res = 700)
rmse_NoLineErrorBar
dev.off()


# save object for combining w/ other figure
# add info about sample size
RMSE_df_med$shapeID <- NA
RMSE_df_med[RMSE_df_med$RMSE_n > 10000, "shapeID"] <- "n ~ 18,000"
RMSE_df_med[RMSE_df_med$RMSE_n < 10000, "shapeID"] <- "n ~ 5,000"

(rmse_NoLineErrorBarNew <- ggplot(data = RMSE_df_med) +
    facet_grid(.~missingness) +
    ggh4x::facet_grid2(
      ~ factor(missingness, levels = c("MCAR", "MNAR"), labels =c("'Missing Completely at Random'", "'Missing Not at Random'")),
      labeller =  label_parsed) +
    geom_linerange(aes(x = propMiss_binned, ymin = quantile_low, ymax = quantile_high, color = modelType), alpha = 1, position = position_dodge(width = .1)) +
    geom_point(aes(x = propMiss_binned, y = RMSE_med, color = modelType, shape = as.factor(shapeID)), alpha = 1, position = position_dodge(width = .1)) +
    #geom_text(aes(x = propMiss_binned, y = 0.8, label = paste0("N=",RMSE_n), group = modelType),  hjust = 0.2, angle =90, position = position_dodge(width = .1) ) + 
    #geom_smooth(aes(x = propMiss_bin, y = RMSE_mean, col = type), method = "lm", se = FALSE) +
    theme_classic() +
    ylab("Root Mean Square Error (RMSE)") +
    scale_shape_manual(values = c(15,0)) +
    xlab("Proportion of Missing Data") + 
    #ylim(c(0,1.25)) + 
    scale_x_continuous(breaks=c(0.2,0.4, 0.6)) +
    
    scale_color_discrete(type = c("#CC79A7", "#E69F00", "#D55E00", "#8c8c8c", "#009E73"),
                         labels = c("Data Augmentation", "Data Deletion-Simple", "Data Deletion-Complete", "Expectation Maximization", "Multiple Imputation")
    )  +
    guides(col = "none", 
           shape = guide_legend(title = "Sample \nSize")#guide_legend(title = "Model Type", position = "right", direction = "vertical", nrow = 5
      ))

saveRDS(rmse_NoLineErrorBarNew, "./figures/RMSEfig_poissSim.rds")

