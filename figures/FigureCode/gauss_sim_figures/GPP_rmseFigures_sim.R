#####
# making figures for GPP forecast RMSE 
#####

# load packages
library(tidyverse)
library(Metrics)
library(magick)
library("grid")
library("ggplotify")

# read in prediction data for au-sable river ------------------------------------------------------------

 #mutate(date = lubridate::as_datetime(date))

# MAR arima data
MCAR_arima <- read.csv("./data/model_results/gauss_sim_randMiss_modelResults_A/AllPreds_arima.csv") %>% 
  rbind(read.csv("./data/model_results/gauss_sim_randMiss_modelResults_B/AllPreds_arima.csv"))
MCAR_arima <- MCAR_arima %>% 
  rename(Estimate = pred, Est.Error = se) %>% 
  mutate(Q2.5 = NA, Q97.5 = NA) %>% 
  relocate(missingprop_autocor, Estimate, Est.Error, Q2.5, Q97.5, date, GPP, missingness, type, sim_no) %>% 
  mutate(date = lubridate::as_datetime(date))

# MAR brms data
MCAR_brms_A <- read.csv("./data/model_results/gauss_sim_randMiss_modelResults_A/brmsResults/AllPreds_brms.csv") #%>% 

MCAR_brms_B <- read.csv("./data/model_results/gauss_sim_randMiss_modelResults_B/brmsResults/AllPreds_brms.csv") 

MCAR_brms <- MCAR_brms_A %>% 
  rbind(MCAR_brms_B) %>% 
  mutate(date = lubridate::as_datetime(date), 
         curSim = NA) %>% 
  rename("sim_no" = "run_no") %>% 
  select(-X)

# MNAR arima data
MNAR_arima <- read.csv("./data/model_results/gauss_sim_minMax_modelResults/AllPreds_arima.csv") 
MNAR_arima <- MNAR_arima %>% 
  rename(Estimate = pred, Est.Error = se) %>% 
  mutate(Q2.5 = NA, Q97.5 = NA) %>% 
  relocate(missingprop_autocor, Estimate, Est.Error, Q2.5, Q97.5, date, GPP, missingness, type, sim_no) %>% 
  mutate(date = lubridate::as_datetime(date),
         missingness = "MNAR")

# MNAR brms data
MNAR_brms <- read.csv("./data/model_results/gauss_sim_minMax_modelResults/AllPreds_brms.csv")%>% 
  mutate(date = lubridate::as_datetime(date),
         missingness = "MNAR", 
         curSim = NA) %>% 
  rename("sim_no" = "run_no")%>% 
  select(-X)

# get real, complete dataset
realData <- MCAR_arima %>% 
  filter(missingprop_autocor == "y_noMiss") %>% 
  mutate(Estimate = NA, 
         Est.Error = NA, 
         missingness = NA, 
         type = "realData", 
         curSim = NA) %>% 
  unique() %>% 
  rename("GPP_actual" = "GPP")


# join all data together
allDat_temp <- rbind(MCAR_arima, MCAR_brms, MNAR_arima, MNAR_brms)

# remove rows w/ no missingness (missingprop_autocor = "y)
allDat_temp <- allDat_temp %>%
  filter(missingprop_autocor != "y_noMiss") %>%
  filter(!is.na(Estimate)) # remove NAs in Estimate column (12/31)

# reformat data -----------------------------------------------------------
# retrieve proportion missing and amount autocor from names
# fix issue with autocorrelation calculation
allDat_temp$missingprop_autocor[str_detect(allDat_temp$missingprop_autocor, pattern = regex("\\...\\.[0-9]$"))] <- 
  str_sub(allDat_temp$missingprop_autocor[str_detect(allDat_temp$missingprop_autocor, pattern = regex("\\...\\.[0-9]$"))], 
          start = 1, end = 30)
allDat_temp$missingprop_autocor[str_detect(allDat_temp$missingprop_autocor, pattern = regex("\\...\\.[0-9]$"))] <- 
  str_sub(allDat_temp$missingprop_autocor[str_detect(allDat_temp$missingprop_autocor, pattern = regex("\\...\\.[0-9]$"))], 
          start = 1, end = 16)

# bin autocorrelation
allDat <- allDat_temp %>% 
  mutate(propMiss = as.numeric(str_split(missingprop_autocor, pattern = "_", simplify = TRUE)[,2]),
         amtAutoCorr = as.numeric(str_split(missingprop_autocor, pattern = "_", simplify = TRUE)[,4]))
allDat[allDat$missingness == "MNAR", "amtAutoCorr"] <- 0
# assign "bins" of autocorrelation
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr <=0.25 & !is.na(allDat$amtAutoCorr), "missingness"] <- "MCAR: low autocorrelation"
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr > 0.25 & allDat$amtAutoCorr < 0.65 &!is.na(allDat$amtAutoCorr), "missingness"] <- "MCAR: medium autocorrelation"
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr  >= 0.65 & !is.na(allDat$amtAutoCorr), "missingness"] <- "MCAR: high autocorrelation"
allDat[allDat$missingness == "MNAR", "missingness"] <- "Missing NOT at Random"
# bin autocorr and propMiss 
allDat$propMiss_bin <- NA
allDat$amtAutoCorr_bin  <- NA

allDat[allDat$propMiss == 0 &  !is.na(allDat$propMiss), "propMiss_bin"] <- 0
allDat[allDat$propMiss > 0 & allDat$propMiss < 0.3 & !is.na(allDat$propMiss), "propMiss_bin"] <- 0.2
allDat[allDat$propMiss >= 0.3 & allDat$propMiss < 0.5 & !is.na(allDat$propMiss), "propMiss_bin"] <- 0.4
allDat[allDat$propMiss >= 0.5 &  !is.na(allDat$propMiss), "propMiss_bin"] <- 0.6

allDat[allDat$amtAutoCorr > 0 & allDat$amtAutoCorr <= 0.25 & !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- "low_autocorr"
allDat[allDat$amtAutoCorr > 0.25 & allDat$amtAutoCorr < 0.65 & !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- "med_autocorr"
allDat[allDat$amtAutoCorr >= 0.65 &  !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <-"high_autocorr"

## fix issue w/ some runs not having real GPP values... can get from other runs (every date should have the same "real" values)
# allDat$GPP <- round(allDat$GPP, 9)
# realGPP <- unique(allDat[,c("date", "GPP", "sim_no")])
# 
# # add back real GPP values to the main data.frame 
# test <- allDat %>% 
#   select(-GPP) %>% 
#   left_join(realGPP, by = c("date", "sim_no"))
# 
# 
rm(MCAR_arima, MCAR_brms, MNAR_arima, MNAR_brms)
gc()


## join the real data to the model predictions
#View(realData)
allDat <- allDat %>% 
 # slice(1:1000) %>%  
  left_join(realData %>% select(date, GPP_actual, sim_no)) 

# # calculate RMSE  ---------------------------------------------------------
RMSE <- allDat %>% 
  group_by(missingprop_autocor, missingness, type, propMiss_bin, amtAutoCorr_bin, sim_no) %>% 
  summarise(RMSE = sqrt(mean((Estimate - GPP)^2, na.rm = TRUE))) 

# add to all data df
allDat <- allDat %>% 
  left_join(RMSE)

# Make figure of predictions ----------------------------------------------
allDat_fig <- allDat %>%
  group_by(date, missingness, type, propMiss_bin) %>%
  summarize(Estimate_mean = mean(Estimate),
            Estimate_sd = sd(Estimate),
            Est.Error_mean = mean(Est.Error)) %>%
  rename(propMiss = "propMiss_bin") %>% 
  mutate(propMiss_fac = factor(propMiss))


RMSE <- RMSE %>% 
  filter(missingness %in% c("Missing NOT at Random", "MCAR: medium autocorrelation")) %>% 
  mutate(missingness = str_replace(missingness, pattern = "MCAR: medium autocorrelation", replacement = "Missing Completely at Random"))

RMSE_errorBar <- RMSE %>% 
  # filter(missingness %in% c("Missing NOT at Random", "MAR: medium autocorrelation"))%>% 
  # mutate(missingness = str_replace(missingness, pattern = "MAR: medium autocorrelation", replacement = "Missing at Random")) %>%
  group_by(missingness, type, propMiss_bin) %>% 
  summarize(RMSE_mean = median(RMSE),
    IQR_high = quantile(RMSE, .75), 
    IQR_low = quantile(RMSE, .25),
    RMSE_n = sum(!is.na(RMSE))) 


## order the missingness type factor
RMSE_errorBar <- RMSE_errorBar %>% 
  mutate(type = factor(type, levels =c("dropNA_simple", "dropNA_complete", "Multiple Imputations",
                       "Kalman Filter", "brms"), ordered = TRUE))
RMSE <- RMSE %>% 
  mutate(type = factor(type, levels =c("dropNA_simple", "dropNA_complete", "Multiple Imputations",
                                       "Kalman Filter", "brms"), ordered = TRUE))


(rmse_NoLineErrorBar <- ggplot(data = RMSE_errorBar) +
    facet_grid(.~missingness) +
    geom_linerange(aes(x = propMiss_bin, ymin = IQR_low, ymax = IQR_high, color = type), alpha = 1, position = position_dodge(width = .1)) +
    geom_point(aes(x = propMiss_bin, y = RMSE_mean, color = type), alpha = 1, position = position_dodge(width = .1)) +
    #geom_smooth(aes(x = propMiss_bin, y = RMSE_mean, col = type), method = "lm", se = FALSE) +
   geom_text(aes(x = propMiss_bin, y = 0.8, label = paste0("N=",RMSE_n), group = type),  hjust = 0.2, angle =90, position = position_dodge(width = .1) ) + 
    theme_classic() +
    ylab("Root Mean Square Error (RMSE)") +
    xlab("Proportion of Missing Data") + 
    #ylim(c(0,1.25)) + 
    scale_x_continuous(breaks=c(0.2,0.4, 0.6)) +
    # scale_color_discrete(type = c("#D55E00","#CC79A7", "#E69F00", "#0072B2","#009E73"),
    #                      labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Kalman Filter", "Multiple Imputation"), 
    # )  +
    scale_color_discrete(type = c( "#E69F00","#D55E00","#009E73", "#0072B2","#CC79A7"),
                         labels = c("Data Deletion-Simple","Data Deletion-Complete",  "Multiple Imputation",  "Kalman Filter","Data Augmentation"), 
    )  +
    guides(col = guide_legend(title = "Model Type", position = "top", direction = "vertical", nrow = 2))
)

png(file = "./figures/RMSE_FullFigure_NoLineWithErrorBar_gaussian_SIM.png", width = 6.5, height = 6, units = "in", res = 700)
rmse_NoLineErrorBar
dev.off()

## save figure to add to empirical RMSE figure
# add sample size info


RMSE_errorBar$shapeID <- NA
RMSE_errorBar[RMSE_errorBar$RMSE_n > 10000, "shapeID"] <- "n ~ 19,000"
RMSE_errorBar[RMSE_errorBar$RMSE_n < 10000, "shapeID"] <- "n ~ 4,000"

(rmse_NoLineErrorBarNew <- ggplot(data = RMSE_errorBar) +
    facet_grid(.~missingness) +
    geom_linerange(aes(x = propMiss_bin, ymin = IQR_low, ymax = IQR_high, color = type), alpha = 1, position = position_dodge(width = .1)) +
    geom_point(aes(x = propMiss_bin, y = RMSE_mean, color = type, shape = as.factor(shapeID)), alpha = 1, position = position_dodge(width = .1)) +
    #geom_smooth(aes(x = propMiss_bin, y = RMSE_mean, col = type), method = "lm", se = FALSE) +
    #geom_text(aes(x = propMiss_bin, y = 0.8, label = paste0("N=",RMSE_n), group = type),  hjust = 0.2, angle =90, position = position_dodge(width = .1) ) + 
    theme_classic() +
    ylab("Root Mean Square Error (RMSE)") +
    xlab("Proportion of Missing Data") + 
    scale_shape_manual(values = c(15,0)) +
    #ylim(c(0,1.25)) + 
    scale_x_continuous(breaks=c(0.2,0.4, 0.6)) +
    # scale_color_discrete(type = c("#D55E00","#CC79A7", "#E69F00", "#0072B2","#009E73"),
    #                      labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Kalman Filter", "Multiple Imputation"), 
    # )  +
    scale_color_discrete(type = c( "#E69F00","#D55E00","#009E73", "#0072B2","#CC79A7"),
                         labels = c("Data Deletion-Simple","Data Deletion-Complete",  "Multiple Imputation",  "Kalman Filter","Data Augmentation"), 
    )  +
    guides(col = "none", 
           shape = guide_legend(title = "Sample \nSize")#guide_legend(title = "Model Type", position = "right", direction = "vertical", nrow = 5
    )
)


saveRDS(rmse_NoLineErrorBarNew, "./figures/RMSEfig_gaussSim.rds")
