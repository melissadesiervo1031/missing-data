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
# read in real, complete dataset
realData <- read.csv("./data/au_sable_river_prepped.csv") %>% 
  mutate(date = lubridate::as_datetime(date))

# MAR arima data
MCAR_arima <- read.csv("./data/model_results/gauss_real_MAR_arima_modResults/au_sable/gauss_auSable_real_MAR_arima_FORECASTpreds.csv") %>% 
  rename(Estimate = pred, Est.Error = se) %>% 
  mutate(Q2.5 = NA, Q97.5 = NA) %>% 
  select(-CurSim, -X, -...1, -...11, -...12, -...13) %>% 
  relocate(missingprop_autocor, Estimate, Est.Error, Q2.5, Q97.5, date, GPP, missingness, type, run_no) %>% 
  mutate(date = lubridate::as_datetime(date))
# MAR brms data
MCAR_brms <- read.csv("./data/model_results/gauss_real_MAR_brms_modResults/auSable/gauss_auSable_real_MAR_brms_FORECASTpreds.csv") %>% 
  mutate(date = lubridate::as_datetime(date)) %>% 
  select(-...1, -X, -...12)

# MNAR arima data
MNAR_arima <- read.csv("./data/model_results/gauss_real_MNAR_arima_modResults/au_sable/gauss_auSable_real_MNAR_arima_FORECASTpreds.csv") %>% 
  rename(Estimate = pred, Est.Error = se) %>% 
  mutate(Q2.5 = NA, Q97.5 = NA) %>% 
  select(-CurSim, -X) %>% 
  relocate(missingprop_autocor, Estimate, Est.Error, Q2.5, Q97.5, date, GPP, missingness, type, run_no) %>% 
  mutate(date = lubridate::as_datetime(date),
         missingness = "MNAR")
# MNAR brms data
MNAR_brms <- read.csv("./data/model_results/gauss_real_MNAR_brms_modResults/auSable/brmspreds.csv")%>% 
  mutate(date = lubridate::as_datetime(date),
         missingness = "MNAR")

# join all data together
allDat_temp <- rbind(MCAR_arima, MCAR_brms, MNAR_arima, MNAR_brms)

# remove rows w/ no missingness (missingprop_autocor = "y)
allDat_temp <- allDat_temp %>% 
  filter(missingprop_autocor != "y") %>% 
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
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr <=0.3 & !is.na(allDat$amtAutoCorr), "missingness"] <- "MCAR: low autocorrelation"
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr > 0.3 & allDat$amtAutoCorr < 0.6 &!is.na(allDat$amtAutoCorr), "missingness"] <- "MCAR: medium autocorrelation"
allDat[allDat$missingness != "MNAR" & allDat$amtAutoCorr  >= 0.6 & !is.na(allDat$amtAutoCorr), "missingness"] <- "MCAR: high autocorrelation"
allDat[allDat$missingness == "MNAR", "missingness"] <- "Missing NOT at Random"
# bin autocorr and propMiss 
allDat$propMiss_bin <- NA
allDat$amtAutoCorr_bin  <- NA

allDat[allDat$propMiss == 0 &  !is.na(allDat$propMiss), "propMiss_bin"] <- 0
allDat[allDat$propMiss > 0 & allDat$propMiss < 0.3 & !is.na(allDat$propMiss), "propMiss_bin"] <- 0.2
allDat[allDat$propMiss >= 0.3 & allDat$propMiss < 0.5 & !is.na(allDat$propMiss), "propMiss_bin"] <- 0.4
allDat[allDat$propMiss >= 0.5 &  !is.na(allDat$propMiss), "propMiss_bin"] <- 0.6

allDat[allDat$amtAutoCorr == 0 &  !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- 0
allDat[allDat$amtAutoCorr > 0 & allDat$amtAutoCorr < 0.375 & !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- 0.25
allDat[allDat$amtAutoCorr >= 0.375 & allDat$amtAutoCorr < 0.625 & !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- 0.5
allDat[allDat$amtAutoCorr >= 0.625 &  !is.na(allDat$amtAutoCorr), "amtAutoCorr_bin"] <- 0.75

## fix issue w/ some runs not having real GPP values... can get from other runs (every date should have the same "real" values)
allDat$GPP <- round(allDat$GPP, 9)
realGPP <- unique(allDat[,c("date", "GPP")])

# add back real GPP values to the main data.frame 

test <- allDat %>% 
  select(-GPP) %>% 
  left_join(realGPP, by = "date")


# calculate RMSE  ---------------------------------------------------------
RMSE <- allDat %>% 
  group_by(missingprop_autocor, missingness, type, propMiss_bin, amtAutoCorr_bin) %>% 
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
    IQR_low = quantile(RMSE, .25)) 


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
    theme_classic() +
    ylab("Root Mean Square Error (RMSE)") +
    xlab("Proportion of Missing Data") + 
    #ylim(c(0,1.25)) + 
    scale_x_continuous(breaks=c(0.2,0.4, 0.6)) +
    scale_color_discrete(type = c("#D55E00","#CC79A7", "#E69F00", "#0072B2","#009E73"),
                         labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Kalman Filter", "Multiple Imp."), 
    )  +
    guides(col = guide_legend(title = "Model Type", position = "top", direction = "vertical", nrow = 2))
)
## get complete time series figure 
aus <- read.csv("./data/au_sable_river_prepped.csv", header=TRUE)

ausNew <- aus %>% 
  mutate(date = as.POSIXct(date)) 

# how many days are in the forecast vs. training data? 
train <- ausNew %>% 
  filter(date < as.POSIXct("2014-01-01T00:00:00Z")) %>% 
  filter(!is.na(GPP))
test <- ausNew %>% 
  filter(date >= as.POSIXct("2014-01-01T00:00:00Z")) %>% 
  filter(!is.na(GPP))

(tsFigGrob <- 
ggplot() + 
  geom_rect(aes(xmin = as.POSIXct("2014-01-01T00:00:00Z"), xmax = as.POSIXct("2014-12-31T00:00:00Z"), 
                ymin = -6, ymax = 3.5), fill = "grey80") +
  geom_line(data = ausNew, aes(x = date, y = GPP)) + 
  geom_rug(data = ausNew[is.na(ausNew$GPP),], aes(date), col = "red") +
  theme_classic() + 
  labs(x = "Year", y = "GPP (scaled)"))

tsPlusRmse_NoLine <- ggpubr::ggarrange(tsFigGrob, rmse_NoLineErrorBar, ncol = 1, nrow = 2, 
                                         heights = c(.5, 1), legend = "bottom", labels = c("A", "B"))



png(file = "./figures/RMSE_FullFigure_NoLineWithErrorBar_gaussian_auSable.png", width = 6.5, height = 8, units = "in", res = 700)
tsPlusRmse_NoLine
dev.off()



# Figure of RMSE variation  -----------------------------------------------
# calculate the width of the IQR for each point shown in the previous figure
RMSE_errorBar <- RMSE_errorBar %>% 
  mutate(IQR_width = IQR_high - IQR_low)

(rmseWidth_fig <- ggplot(data = RMSE_errorBar) +
    facet_grid(.~missingness) +
    geom_linerange(aes(x = propMiss_bin, ymin = 0, ymax = IQR_width, color = type), alpha = 1, position = position_dodge(width = .1), lwd = 2) +
    theme_classic() +
    ylab("Width of RMSE Inter-Quartile Range") +
    xlab("Proportion of Missing Data") + 
    #ylim(c(0,1.25)) + 
    scale_x_continuous(breaks=c(0.2,0.4, 0.6)) +
    scale_color_discrete(type = c("#D55E00","#CC79A7", "#E69F00", "#0072B2","#009E73"),
                         labels = c("Data Deletion-Complete", "Data Augmentation", "Data Deletion-Simple", "Kalman Filter", "Multiple Imp."), 
    )  +
    guides(col = guide_legend(title = "Model Type", position = "top", direction = "vertical", nrow = 2))
)


png(file = "./figures/RMSE_IQR_width_gaussian_auSable.png", width = 6, height = 5, units = "in", res = 700)
rmseWidth_fig
dev.off()

